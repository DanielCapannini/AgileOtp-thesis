//-----------------------------------------------------------------------------
//  File: Knapsack.h                                                       
//                                                                           
//  Project: Knapsack solver based on a trivial Dynamic Programing Procedure                    
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 07.10.2011                                                  
//-----------------------------------------------------------------------------
//

#define  Infkp  1000000000.
#define  ReMap  1

class Knapsack
{
private:

public:

	// Input
	long n;      // Numero di oggetti
	long nf;     // Numero di oggetti fixed to 1
	long W;      // Capacita' del knapsack
	long *w;     // Vettore dei pesi degli oggetti
	long *wf;    // Vettore dei pesi degli oggetti fixed to 1 
	double *p;   // Vettore dei profitti (degli oggetti)
	long *ind;   // Indice oggetto (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)
	long *indf;  // Indice oggetto fixed to 1 (un appoggio per l'utente per mappare chi è stato inserito nel knapsack)
	long fsolver; // fsolver=0 use Dynamic Programing (default); fsolver=1 use Cplex
	long fmalloc; // fmalloc=0 data structure not allocated; fmalloc=1 allocated

	// Output
	double zopt;  // Valore della soluzione ottima
	double zf;    // Valore della soluzione dovuta agli oggetti fixed to 1
	long *x;      // Vettore soluzione: x[i]=0 se non e' in soluzione; x[i]=1 altrimenti
	long *fcover; // Vettore cover: fcover[i]=0 se non e' nella cover; fcover[i]=1 altrimenti

	// Auxiliary
	long maxN;
	long maxW;
	double **ztab;
	long **ftab;
	long *map;
	long *mapw;

	// Constructor and Destructor
	Knapsack(void);
	~Knapsack(void);

	// Malloc: allocate the required data structure                                    
	int Knapsack::Malloc(long maxn, long maxw);

	// Solve the knapsack problem loaded
	int Knapsack::Solve(void);
	int Knapsack::SolveReMap(void);

};
