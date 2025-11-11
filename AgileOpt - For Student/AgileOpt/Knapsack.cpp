//-----------------------------------------------------------------------------
//  File: Knapsack.cpp                                                       
//                                                                           
//  Project: Knapsack solver based on a trivial Dynamic Programing Procedure                    
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 07.10.2011                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>

#include "knapsack.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
Knapsack::Knapsack(void)
{
	// Input
	n = 0;     // Numero di oggetti
	nf = 0;    // Numero di oggetti fixed to 1
	W = 0;     // Capacita' del knapsack
	w = NULL;  // Vettore dei pesi degli oggetti
	wf = NULL; // Vettore dei pesi degli oggetti fixed to 1
	p = NULL;  // Vettore dei profitti (degli oggetti)
	ind = NULL;  // Indice oggetto (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)
	indf = NULL; // Indice oggetto fixed to 1 (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)
	fsolver = -1; // Solver=0 use Dynamic Programing (default); Solver=1 use Cplex
	fmalloc = 0;  // fmalloc=0 data structure not allocated; fmalloc=1 allocated

	// Output
	zopt = 0;  // Valore della soluzione ottima
	zf = 0;    // Valore della soluzione relativa agli oggetti fixed to 1
	x = NULL;  // Vettore soluzione: x[i]=0 se non e' in soluzione; x[i]=1 altrimenti
	fcover = NULL;  // Vettore cover: fcover[i]=0 se non e' nella cover; fcover[i]=1 altrimenti

	// Auxiliary
	ztab = NULL;
	ftab = NULL;
	map = NULL;
	mapw = NULL;


}


//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
Knapsack::~Knapsack(void)
{
	long i;

	// Input
	n = 0;    // Numero di oggetti
	nf = 0;   // Numero di oggetti fixed to 1
	W = 0;    // Capacita' del knapsack
	if (w!=NULL) delete [] w; // Vettore dei pesi degli oggetti
	if (wf!=NULL) delete [] wf; // Vettore dei pesi degli oggetti fixed to 1
	if (p!=NULL) delete [] p; // Vettore dei profitti (degli oggetti)
	if (ind!=NULL) delete [] ind; // Indice oggetto (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)
	if (indf!=NULL) delete [] indf; // Indice oggetto fixed to 1 (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)

	// Output
	zopt = 0; // Valore della soluzione ottima
	zf = 0;   // Valore della soluzione relativa agli oggetti fixed to 1
	if (x!=NULL) delete [] x; // Vettore soluzione: x[i]=0 se non e' in soluzione; x[i]=1 altrimenti
	if (fcover!=NULL) delete [] fcover; // Vettore cover: fcover[i]=0 se non e' nella cover; fcover[i]=1 altrimenti

	// Auxiliary
	if (fmalloc>0)
	{
		for (i=0; i<=maxN; i++)
		{
			delete [] ztab[i];
			delete [] ftab[i];
		}
		delete [] ztab;
		delete [] ftab;
		delete [] map;
		delete [] mapw;
	}

}


//-----------------------------------------------------------------------------
//  Malloc: allocate the required data structure                                    
//-----------------------------------------------------------------------------
//
int Knapsack::Malloc(long maxn, long maxw)
{
	long i;

	maxN = maxn;
	maxW = maxw;

	w = new long [maxn];
	wf = new long [maxn];
	p = new double [maxn];
	x = new long [maxn];
	fcover = new long [maxn];
	ind = new long [maxn];
	indf = new long [maxn];

	ztab = new double* [maxn+1];
	ftab = new long* [maxn+1];
	for (i=0; i<=maxn; i++)
	{
		ztab[i] = new double [maxw+1];
		ftab[i] = new long [maxw+1];
	}
	map = new long [maxw+1];
	mapw = new long [maxw+1];

	fmalloc=1;

	return 0;
}


//-----------------------------------------------------------------------------
//  Solver: Dynamic Programming Procedure                                    
//-----------------------------------------------------------------------------
//
int Knapsack::Solve(void)
{
	long i;
	long ww;
	double aux;

	// Trivial case
	if (n==0)
		return 0;
	if (n==1)
	{
		if (w[0]>W)
		{
			x[0]=0;
			zopt=0.0;
		}
		else
		{
			x[0]=1;
			zopt=p[0];
		}
		return 0;
	}

	if (fmalloc==0)
	{
		// Allocazione della "tabella" stadi/stati
		ztab = new double* [n+1];
		ftab = new long* [n+1];
#pragma omp parallel for 
		for (i=0; i<=n; i++)
		{
			ztab[i] = new double [W+1];
			ftab[i] = new long [W+1];
		}
	}

	// Passo base
#pragma omp parallel for 
	for (ww=0; ww<=W; ww++)
	{
		ztab[0][ww]=0.0;
		ftab[0][ww]=0;
	}

	// Sviluppo stadi da i=1 a n
	for (i=0; i<n; i++)
#pragma omp parallel for private(aux)
		for (ww=0; ww<=W; ww++)
		{
			if (ww<w[i])
			{
				ztab[i+1][ww]=ztab[i][ww];
				ftab[i+1][ww]=0;
			}
			else
			{
				aux=ztab[i][ww-w[i]]+p[i];
				if (aux>ztab[i][ww])
				{
					ztab[i+1][ww]=aux;
					ftab[i+1][ww]=1;
				}
				else
				{
					ztab[i+1][ww]=ztab[i][ww];
					ftab[i+1][ww]=0;
				}
			}
		}

	//for (ww=0;ww<=W;ww++)
	//	printf("ztab[n][%d]=%lf\n",ww,ztab[n][ww]);

	zopt=ztab[n][W];
	ww=W;
	for (i=n;i>0;i--)
	{
		if (ftab[i][ww]>0)
		{
			x[i-1]=1;
			//printf("  item: %d\n",i-1);
			ww=ww-w[i-1];
		}
		else
		{
			x[i-1]=0;
		}
	}

	if (fmalloc==0)
	{
		// Deallocazione della "tabella" stadi/stati
#pragma omp parallel for 
		for (i=0;i<=n;i++)
		{
			delete [] ztab[i];
			delete [] ftab[i];
		}
		delete [] ztab;
		delete [] ftab;
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  Solver: Dynamic Programming Procedure with ReMapping 
//          (versione simile a quella a liste di Bellman)
//-----------------------------------------------------------------------------
//
int Knapsack::SolveReMap(void)
{
	long i;
	long w1,w2,ww;
	double aux;
	long nmap,nmap1;

	// Trivial case
	if (n==0)
		return 0;
	if (n==1)
	{
		if (w[0]>W)
		{
			x[0]=0;
			zopt=0.0;
		}
		else
		{
			x[0]=1;
			zopt=p[0];
		}
		return 0;
	}

	// Sum of the item smaller than W
	ww=0;
	for (i=0; i<n; i++)
		ww += w[i];

	if (ww<=W)
	{
		zopt = 0.0;
		for (i=0; i<n; i++)
		{
			zopt += p[i];
			x[i] = 1;
		}
		return 0;
	}

	// Allocate data structure if required
	if (fmalloc==0)
	{
		// Allocazione della "tabella" stadi/stati
		ztab = new double* [n+1];
		ftab = new long* [n+1];
//#pragma omp parallel for 
		for (i=0; i<=n; i++)
		{
			ztab[i] = new double [W+1];
			ftab[i] = new long [W+1];
		}

		nmap = 0;
		map = new long [W+1];
		mapw = new long [W+1];
	}

	// Optimize by dynamic programming 
//#pragma omp parallel for 
	for (ww=0; ww<=W; ww++)
	{
		map[ww] = -1;
		mapw[ww] = -1;
	}

	// Passo base
//#pragma omp parallel for 
//	for (ww=0; ww<=W; ww++)
//	{
//		ztab[0][ww]=0.0;
//		ftab[0][ww]=0;
//	}
	map[0] = 0;  // map[w]: il peso w=0 lo trovi in posizione map[w]=0
	mapw[0] = 0; // mapw[i]: in posizione i=0 trovi il peso mapw[i]=0
	nmap = 1;    // UN peso è stato inserito nella lista
	ztab[0][0] = 0.0;
	ftab[0][0] = 0;

	// Sviluppo stadi da i=1 a n
	for (i=0; i<n; i++)
	{
//#pragma omp parallel for private(aux)
		for (ww=0; ww<nmap; ww++)
		{
			ztab[i+1][ww]=ztab[i][ww];
			ftab[i+1][ww]=0;
		}
		nmap1=nmap;
		for (ww=0; ww<nmap1; ww++)
		{
			w1 = mapw[ww] + w[i];
			if (w1>W) 
				continue;

			if (map[w1]<0)
			{
				mapw[nmap] = w1;
				map[w1] = nmap;
				ztab[i+1][nmap]=ztab[i][ww]+p[i];
				ftab[i+1][nmap]=1;
				nmap++;
			}
			else
			{
				w2 = map[w1];
				aux=ztab[i][ww]+p[i];
				if (aux>ztab[i][w2])
				{
					ztab[i+1][w2]=aux;
					ftab[i+1][w2]=1;
				}
			}
		}
	}

	//for (ww=0;ww<=W;ww++)
	//	printf("ztab[n][%d]=%lf\n",ww,ztab[n][ww]);

	zopt = -Infkp;
	for (ww=0; ww<nmap; ww++)
	{
		if (ztab[n][ww]>zopt)
		{
			zopt = ztab[n][ww];
			w1 = ww;
		}
	}
	for (i=n;i>0;i--)
	{
		if (ftab[i][w1]>0)
		{
			x[i-1] = 1;
			//printf("  item: %d\n",i-1);
			ww = mapw[w1]-w[i-1];
			w1 = map[ww];
		}
		else
		{
			x[i-1]=0;
		}
	}

	// Deallocazione della "tabella" stadi/stati
	if (fmalloc==0)
	{
//#pragma omp parallel for 
		for (i=0;i<=n;i++)
		{
			delete [] ztab[i];
			delete [] ftab[i];
		}
		delete [] ztab;
		delete [] ftab;
		delete [] map;
		delete [] mapw;
	}

	return 0;
}
