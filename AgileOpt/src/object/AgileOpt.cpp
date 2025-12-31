//-----------------------------------------------------------------------------
//  File: AgileOpt.cpp                                                       
//                                                                           
//  Project: Optimize the Agile Project Management Schedule                  
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 22.08.2011                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../lib/AgileOpt.h"
#include "../lib/CplexObj.h"
//#include "Knapsack.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
AgileOpt::AgileOpt(void)
{
	ferr = fopen("Agile.err","w");

	n = 0;       // Number of User Stories
	m = 0;       // Number of Sprints
	a = NULL;    // a[j] is the affinity of story j
	u = NULL;    // u[j] is the utility of story j
	ur = NULL;   // u[j]*rcr[j] is the "base" weighted utility of story j
	p = NULL;    // p[j] is the number of story points of story j
	pr = NULL;   // p[j]*run[j] is the weight of story j
	pmax = NULL; // pmax[i] is the capacity of sprint i, measured in story points
	rcr = NULL;  // rcr[j] is the criticality risk of story j
	run = NULL;  // run[j] is the uncertainty risk of story j
	nY = NULL;   // nY[j] is the number of stories similar to story j
	Y = NULL;    // Y[j][k] is the k-th story similar to story j
	nD = NULL;   // nD[j] is the number of stories the story j depends on
	D = NULL;    // D[j][k] is the k-th story the story j depends on
	nUOR = NULL; // nUOR[j] is the number of stories having an OR dependency with story j
	UOR = NULL;  // UOR[j][k] is the k-th story having an OR dependency with story j
	nUAND = NULL;// nUAND[j] is the number of stories having an AND dependency with story j
	UAND = NULL; // UAND[j][k] is the k-th story having an AND dependency with story j
	FDepA = NULL;// FDepA[j] is 1 if the story j is involved in some dependency (Active)
	FDepP = NULL;// FDepP[j] is 1 if the story j is involved in some dependency (Passive)

	cfg.TimeLimit = 10.0*24.0*3600.0;  // 10 gg di runt-time
	cfg.Sentinel = 0;
	cfg.Pred = 0;
	cfg.Cover = 0;
	cfg.Lifting = 0;
	cfg.KnapSol = 0;
	cfg.AddLB = 0;
	cfg.MaxCuts = 100000;
	Zopt = 0.0;
	Gap = 0.0;
	Nodes = 0;
	Cuts = 0;
	NSprintsUsed = 0;
	PercentageUtilization = 0;
	DeviationRisk = 0.0;
	DeviationUncertaintyRisk = 0.0;
	HalfUtilitiSprint = -1;
}

long CountOccupiedSprints(const CplexObj& LP, long n, long m)
{
    long i, j;
    long count = 0;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (LP.xt[i*n + j] > 0.5) {
                count++;
                break;
            }
        }
    }
    return count;
}

double AverageSprintUtilization(const CplexObj& LP, long n, long m, const double* p, const double* run, const double* pmax)
{
    long i, j;
    long occupied = 0;
    double avg = 0.0;

    for (i = 0; i < m; i++) {
        double used_points = 0.0;
        bool used = false;

        for (j = 0; j < n; j++) {
            if (LP.xt[i*n + j] > 0.5) {
                used = true;
                used_points += p[j] * run[j];
            }
        }

        if (used) {
            occupied++;
            avg += used_points / pmax[i];
        }
    }

    return (occupied > 0) ? avg / occupied : 0.0;
}

double SprintRiskStdDev(const CplexObj& LP, long n, long m, const double* rcr)
{
    long i, j;

    double *R = new double[m];
    for (i = 0; i < m; i++)
        R[i] = 0.0;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (LP.xt[i*n + j] > 0.5)
                R[i] += rcr[j];

    double mean = 0.0;
    for (i = 0; i < m; i++)
        mean += R[i];
    mean /= m;

    double var = 0.0;
    for (i = 0; i < m; i++)
        var += (R[i] - mean) * (R[i] - mean);

    delete [] R;

    return sqrt(var / m);
}

double SprintUncertaintyStdDev(const CplexObj& LP, long n, long m, const double* run)
{
    long i, j;

    double *U = new double[m];
    for (i = 0; i < m; i++)
        U[i] = 0.0;

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (LP.xt[i*n + j] > 0.5)
                U[i] += run[j];

    double mean = 0.0;
    for (i = 0; i < m; i++)
        mean += U[i];
    mean /= m;

    double var = 0.0;
    for (i = 0; i < m; i++)
        var += (U[i] - mean) * (U[i] - mean);

    delete [] U;

    return sqrt(var / m);
}

long SprintHalfUtility(const CplexObj& LP, long n, long m, const double* u) 
{
    long i, j;

    double U_tot = 0.0;
    for (j = 0; j < n; j++)
        U_tot += u[j];

    double U_cum = 0.0;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (LP.xt[i*n + j] > 0.5)
                U_cum += u[j];
        }

        if (U_cum >= U_tot / 2.0)
            return i;
    }

    return -1;  // mai raggiunta
}


//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
AgileOpt::~AgileOpt(void)
{
	long j;

	if (a!=NULL) {delete [] a; a = NULL;}    
	if (u!=NULL) {delete [] u; u = NULL;}
	if (ur!=NULL) {delete [] ur; ur = NULL;}
	if (p!=NULL) {delete [] p; p = NULL;}
	if (pr!=NULL) {delete [] pr; pr = NULL;}
	if (pmax!=NULL) {delete [] pmax; pmax = NULL;}
	if (rcr!=NULL) {delete [] rcr; rcr = NULL;}
	if (run!=NULL) {delete [] run; run = NULL;}

	if (nY!=NULL)
	{
		for (j=0; j<n; j++) 
			if (Y[j]!=NULL) {delete [] Y[j]; Y[j] = NULL;}
		delete [] Y;
		delete [] nY;
	}
	if (nD!=NULL)
	{
		for (j=0; j<n; j++) 
			if (D[j]!=NULL) {delete [] D[j]; D[j] = NULL;}
		delete [] D;
		delete [] nD;
	}
	if (nUOR!=NULL)
	{
		for (j=0; j<n; j++) 
			if (UOR[j]!=NULL) {delete [] UOR[j]; UOR[j] = NULL;}
		delete [] UOR;
		delete [] nUOR;
	}

	if (nUAND!=NULL)
	{
		for (j=0; j<n; j++) 
			if (UAND[j]!=NULL) {delete [] UAND[j]; UAND[j] = NULL;}
		delete [] UAND;
		delete [] nUAND;
	}

	if (FDepA!=NULL) {delete [] FDepA; FDepA = NULL;}
	if (FDepP!=NULL) {delete [] FDepP; FDepP = NULL;}

	n = 0;       
	m = 0;       
	
	fclose(ferr);

}

//-----------------------------------------------------------------------------
//  Read the instance from Text File "name"                                    
//-----------------------------------------------------------------------------
//
int AgileOpt::ReadData(const char *name)
{
	FILE *fin;
	long nr,nr1;
	long i,j,i1,j1,i2;
	long nlist;
	long *list;
	double aux;
	char Dep[8];

	fin=fopen(name,"r");

	// Read the number of user stories and of sprints
	fscanf(fin,"%d %d",&n,&m);

	original_m = m;

	// Read Utility
	u = new double [n];
	for (j=0; j<n; j++)
		fscanf(fin,"%lf%",&(u[j]));

	// Read Criticality Risk
	rcr = new double [n];
	ur = new double [n];
	for (j=0; j<n; j++)
	{
		fscanf(fin,"%lf%",&(rcr[j]));
		ur[j] = u[j]*rcr[j];
	}
	// Read Uncertainty Risk
	run = new double [n];
	for (j=0; j<n; j++)
		fscanf(fin,"%lf%",&(run[j]));

	// Read Number of Story Points
	p = new double [n];
	pr = new double [n];
	for (j=0; j<n; j++)
	{
		fscanf(fin,"%lf%",&(p[j]));
		pr[j]=p[j]*run[j];
	}

	// Read Affinity
	a = new double [n];
	nY = new long [n];
	Y = new long* [n];
	list = new long [n];
	for (j=0; j<n; j++)
	{
		a[j] = 0.0;
		nY[j] = 1;
		Y[j] = new long [n];
		Y[j][0] = j;  // We suppose that Yj containts at least j...
	}
	fscanf(fin,"%d",&nr);
	for (i=0; i<nr; i++)
	{
		nlist=0;
		//fscanf(fin,"%lf %d",&aux,&nr1);
		fscanf(fin,"%lf",&aux);
		GetList(fin,&nlist,list);
		//for (i1=0; i1<nr1; i1++)
		//{
		//	fscanf(fin,"%d",&j1);
		//	a[j1]=aux;
		//	list[nlist]=j1;
		//	nlist++;
		//}
		for (i1=0; i1<nlist; i1++)
		{
			j1=list[i1];
			nY[j1]=nlist;
			a[j1]=aux;
			for (i2=0; i2<nlist; i2++)
				Y[j1][i2]=list[i2];
		}
	}
	delete [] list;

	// Read AND/OR Dependancy
	nD = new long [n];
	D = new long* [n];
	nUOR = new long [n];
	UOR = new long* [n];
	nUAND = new long [n];
	UAND = new long* [n];
	FDepA = new long [n];
	FDepP = new long [n];
	for (j=0; j<n; j++)
	{
		nD[j] = 0;
		D[j] = new long [n];
		nUOR[j] = 0;
		UOR[j] = new long [n];
		nUAND[j] = 0;
		UAND[j] = new long [n];
		FDepA[j] = 0;
		FDepP[j] = 0;
	}
	fscanf(fin,"%d",&nr);
	for (i=0; i<nr; i++)
	{
		fscanf(fin,"%d %s",&j1,Dep);
		GetList(fin,&(nD[j1]),D[j1]);

		FDepA[j1] = 1;
		if (strcmp(Dep,"OR")==0)
		{
			nUOR[j1] = nD[j1];
			for (i1=0; i1<nD[j1]; i1++)
			{
				UOR[j1][i1] = D[j1][i1];
				FDepP[D[j1][i1]] = 1;
			}
		}
		else if (strcmp(Dep,"AND")==0)
		{
			nUAND[j1] = nD[j1];
			for (i1=0; i1<nD[j1]; i1++)
			{
				UAND[j1][i1] = D[j1][i1];
				FDepP[D[j1][i1]] = 1;
			}
		}
		else 
		{
			printf("\nERROR... while reading dependency\n");
			fprintf(ferr,"\nERROR... while reading dependency\n");
			exit(1);
		}
	}

	// Read Sprint Capacities
	//m=10;
	pmax = new double [m];
	for (i=0; i<m; i++)
		fscanf(fin,"%lf",&(pmax[i]));

 	fclose(fin);

	return 0;
}

void SkipList(FILE *fin)
{
    long nlist;
    fscanf(fin, "%d", &nlist);
    long dummy;
    for (long k = 0; k < nlist; k++)
        fscanf(fin, "%ld", &dummy);
}


//-----------------------------------------------------------------------------
//  Read the instance from Text File "name"                                    
//-----------------------------------------------------------------------------
//
int AgileOpt::ReadDataWithoutConstraints(const char *name)
{
	FILE *fin;
	long nr,nr1;
	long i,j,i1,j1,i2;
	long nlist;
	long *list;
	double aux;
	char Dep[8];

	fin=fopen(name,"r");

	// Read the number of user stories and of sprints
	fscanf(fin,"%d %d",&n,&m);

	original_m = m;

	// Read Utility
	u = new double [n];
	for (j=0; j<n; j++)
		fscanf(fin,"%lf%",&(u[j]));

	// Read Criticality Risk
	rcr = new double [n];
	ur = new double [n];
	for (j=0; j<n; j++)
	{
		fscanf(fin,"%lf%",&(rcr[j]));
		ur[j] = u[j]*rcr[j];
	}
	// Read Uncertainty Risk
	run = new double [n];
	for (j=0; j<n; j++)
		fscanf(fin,"%lf%",&(run[j]));

	// Read Number of Story Points
	p = new double [n];
	pr = new double [n];
	for (j=0; j<n; j++)
	{
		fscanf(fin,"%lf%",&(p[j]));
		pr[j]=p[j]*run[j];
	}

	// Read Affinity
	a = new double [n];
	nY = new long [n];
	Y = new long* [n];
	list = new long [n];
	for (j=0; j<n; j++)
	{
		a[j] = 0.0;
		nY[j] = 1;
		Y[j] = new long [n];
		Y[j][0] = j;  // We suppose that Yj containts at least j...
	}
	fscanf(fin,"%d",&nr);
	for (i=0; i<nr; i++)
	{
		nlist=0;
		//fscanf(fin,"%lf %d",&aux,&nr1);
		fscanf(fin,"%lf",&aux);
		GetList(fin,&nlist,list);
		//for (i1=0; i1<nr1; i1++)
		//{
		//	fscanf(fin,"%d",&j1);
		//	a[j1]=aux;
		//	list[nlist]=j1;
		//	nlist++;
		//}
		for (i1=0; i1<nlist; i1++)
		{
			j1=list[i1];
			nY[j1]=nlist;
			a[j1]=aux;
			for (i2=0; i2<nlist; i2++)
				Y[j1][i2]=list[i2];
		}
	}
	delete [] list;

	// Read AND/OR Dependancy
	nD = new long [n];
	D = new long* [n];
	nUOR = new long [n];
	UOR = new long* [n];
	nUAND = new long [n];
	UAND = new long* [n];
	FDepA = new long [n];
	FDepP = new long [n];
	for (j=0; j<n; j++)
	{
		nD[j] = 0;
		D[j] = new long [n];
		nUOR[j] = 0;
		UOR[j] = new long [n];
		nUAND[j] = 0;
		UAND[j] = new long [n];
		FDepA[j] = 0;
		FDepP[j] = 0;
	}
	fscanf(fin,"%d",&nr);
	for (i=0; i<nr; i++)
	{
		fscanf(fin,"%d %s",&j1,Dep);
		GetList(fin,&(nD[j1]),D[j1]);

		/*
		FDepA[j1] = 1;
		if (strcmp(Dep,"OR")==0)
		{
			nUOR[j1] = nD[j1];
			for (i1=0; i1<nD[j1]; i1++)
			{
				UOR[j1][i1] = D[j1][i1];
				FDepP[D[j1][i1]] = 1;
			}
		}
		else if (strcmp(Dep,"AND")==0)
		{
			nUAND[j1] = nD[j1];
			for (i1=0; i1<nD[j1]; i1++)
			{
				UAND[j1][i1] = D[j1][i1];
				FDepP[D[j1][i1]] = 1;
			}
		}
		else 
		{
			printf("\nERROR... while reading dependency\n");
			fprintf(ferr,"\nERROR... while reading dependency\n");
			exit(1);
		}*/
	}

	// Read Sprint Capacities
	//m=10;
	pmax = new double [m];
	for (i=0; i<m; i++)
		fscanf(fin,"%lf",&(pmax[i]));

 	fclose(fin);

	return 0;
}


//-----------------------------------------------------------------------------
//  Utility for reading a list of integer values of unknown lenght                                    
//-----------------------------------------------------------------------------
//
long AgileOpt::GetList(FILE *fin, long *nlist, long *list)
{
	static char line[255];
	char cc;
	long k,h,status;
	int bol;

	h=0;
loop0:

	k=0;
	bol=0;
loop1:
	cc=fgetc(fin);
	if ((cc>='0')&&(cc<='9'))
	{
		bol=1;
		line[k]=cc;
		k++;
		goto loop1;
	}

	if (bol<1)
	{
		if (cc!='\n')
			goto loop0;
	}
	else
	{
		line[k]='\0';
		list[h]=atoi(line);
		h++;
		if (cc!='\n')
			goto loop0;
	}

	*nlist=h;

	return k;
}

//-----------------------------------------------------------------------------
//  Optimize the Agile Schedule                                    
//-----------------------------------------------------------------------------
//
int AgileOpt::Optimize(void)
{
	long i,j,k,h,j1,kk;
	long kount;
	double obj;
	double ptot;
	double objtot;
	double objval;
	long nYMap;
	long *YMap;

	// Reductions
	// Reduction();

	// If required compute a heuristic solution
	if (cfg.AddLB)
		OptimizeHeu();

	// We allocate Cplex here
	CplexObj LP;

	// Auxiliary array to map variables "yij"
	YMap = new long [n*m];
	k = 0;
	nYMap = 0;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		{
			if ((nY[j]>1)&&(a[j]>Prec))
			//if (a[j]>-Inf)
			{
				YMap[k]=nYMap;
				nYMap++;
			}
			else
				YMap[k]=-1;
			k++;
		}

	// Set problem size
	//LP.ncols = 2*n*m;
	//LP.nrows = m + n + 2*n*m;  // OR/AND constraints are added later
	//LP.nz = (n*m)*(3+n) + (n*m)*(1+n);  // We overestimate |Yj|=n

	LP.ncols = n*m+nYMap;
	LP.nrows = m + n + 2*nYMap;  // OR/AND constraints are added later
	LP.nz = (n*m)*(3+n) + (nYMap)*(1+n);  // We overestimate |Yj|=n
	if (cfg.Sentinel)
	{
		LP.ncols += m;
		LP.nrows += m;
		LP.nz += n+m;
	}
	LP.nsos = n;
	LP.nsosnz = n*m;
	LP.MallocCols();

	// Load Columns/Variables xij in the Cplex data structure
	k = 0;
	kount = 0;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		{
			//k = i*n+j;
			LP.obj[k] = (original_m-i)*u[j]*rcr[j];
			//LP.obj[k] = (i+1)*u[j]*rcr[j];

			LP.matbeg[k] = kount;

			LP.matind[kount] = i;
			// LP.matval[kount] = p[j]*run[j]; // or pr[j]
			LP.matval[kount] = pr[j];
			kount++;

			LP.matind[kount] = m+j;
			LP.matval[kount] = 1.;
			kount++;

			for (h=0; h<nY[j]; h++)
			{
				j1 = Y[j][h];
				if ((YMap[i*n+j1]>=0)&&(j1!=j))
				{
					//LP.matind[kount] = m+n+i*n+j1;
					LP.matind[kount] = m+n+YMap[i*n+j1];
					LP.matval[kount] = -1.;
					kount++;
				}
			}

			if (YMap[i*n+j]>=0)
			{
				//LP.matind[kount] = m+n+n*m+i*n+j;
				LP.matind[kount] = m+n+nYMap+YMap[i*n+j];
				LP.matval[kount] = -(double)(nY[j]-1);
				//LP.matval[kount] = -(double)(nY[j]);
				kount++;
			}

			// Sentinel Constraint
			if (cfg.Sentinel)
			{
				LP.matind[kount] = m + n + 2*nYMap + i;
				LP.matval[kount] = 1.0;
				kount++;
			}

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;  // prima 1
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'B';

			LP.lb[k] = 0.0;
			LP.ub[k] = 1.0;
			k++;
		}

	// Load Columns/Variables yij in the Cplex data structure
	kk = 0;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
		{
			if (YMap[kk]<0)
			{
				kk++;
				continue;
			}

			//k = i*n+j;
			if (nY[j]<2)
				LP.obj[k] = 0.0;
			else
				LP.obj[k] = (original_m-i)*u[j]*a[j]/(nY[j]-1);

			//LP.obj[k] = (i+1)*a[j]/nY[j];
			//LP.obj[k] = 0.0;

			LP.matbeg[k] = kount;

			//LP.matind[kount] = m+n+i*n+j;
			LP.matind[kount] = m+n+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			//LP.matind[kount] = m+n+n*m+i*n+j;
			LP.matind[kount] = m+n+nYMap+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'C';
			//LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

	if (cfg.Sentinel)
	{
		// Load Columns/Variables related to Sentinel constraints in the Cplex data structure
		kk = 0;
		for (i=0; i<m; i++)
		{
			LP.obj[k] = 0.0;

			LP.matbeg[k] = kount;

			LP.matind[kount] = m + n + 2*nYMap + i;
			LP.matval[kount] = -1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 1000; // Prima 0
			//LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}
	}

	LP.matbeg[k]=kount;

	// Load Constraints in the Cplex data structure
	kount=0;
	for (i=0; i<m; i++)
	{
		LP.rhs[kount]=pmax[i];
		LP.sense[kount]='L';
		kount++;
	}

	for (j=0; j<n; j++)
	{
		LP.rhs[kount]=1.0;
		LP.sense[kount]='E';
		kount++;
	}

	for (k=0; k<nYMap; k++)
	{
		LP.rhs[kount]=0.0;
		LP.sense[kount]='L';
		kount++;
	}

	for (k=0; k<nYMap; k++)
	{
		LP.rhs[kount]=0.0;
		LP.sense[kount]='L';
		kount++;
	}

	if (cfg.Sentinel)
	{
		for (i=0; i<m; i++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='E';
			kount++;
		}
	}

	//for (i=0; i<m; i++)
	//	for (j=0; j<n; j++)
	//	{
	//		LP.rhs[kount]=0.0;
	//		LP.sense[kount]='L';
	//		kount++;
	//	}

	//for (i=0; i<m; i++)
	//	for (j=0; j<n; j++)
	//	{
	//		LP.rhs[kount]=0.0;
	//		LP.sense[kount]='L';
	//		kount++;
	//	}

    // Define SOS
	k = 0;
	for (j=0; j<n; j++)
	{
		LP.typesos[j] = CPX_TYPE_SOS1;
		LP.sosbeg[j] = k;
		for (i=0; i<m; i++)
		{
			LP.sosind[k] = i*n+j;
			// LP.soswt[k] = (double)(m-i);
			LP.soswt[k] = (double)(i);
			k++;
		}
	}
	LP.sosbeg[n]=k;

	// Maximization problem
	LP.minmax=-1;

	// Load Problem
	LP.CopyLP();

	// Add Integer Sentinel Constraints
	//LP.Sentinel(n,m);

	// Set MIP
	LP.SetMIP(cfg.TimeLimit);

	// Setup Custom Cutting Plane
	LP.CuttingPlane(n,m,ur,pr,pmax,nY,Y,nUOR,UOR,nUAND,UAND,FDepA,FDepP,cfg.Pred,cfg.Cover,cfg.Lifting,cfg.KnapSol,cfg.MaxCuts,&Cuts);

	if (cfg.AddLB>0)
		LP.SetLB(Zheu);

	// Optimize
	LP.SolveMIP(&Zopt,&Gap,&Nodes,&Cuts);

	// Read the solution
	printf("\n\nSoluzione: \n\n");
	objval=0.0;
	k=0;
	kk=n*m;

	NSprintsUsed = 0;
	PercentageUtilization = 0;
	DeviationRisk = 0.0;
	DeviationUncertaintyRisk = 0.0;
	HalfUtilitiSprint = -1;
	double *sprint_risks = new double[m];
	double *sprint_uncertainties = new double[m];
	for (i = 0; i < m; i++) {
		sprint_risks[i] = 0.0;
		sprint_uncertainties[i] = 0.0;
	}

	double halftotal_utility = 0.0;
	for (j=0; j<n; j++)
		halftotal_utility += u[j];
	halftotal_utility /= 2.0;

	int count = 0;
	for (i=0; i<m; i++)
	{
		ptot=0.0;
		objtot=0.0;
		printf("Sprint %d (pmax=%lf): ",i,pmax[i]);
		double sprint_utilization = 0.0;

		double avg_sprint_risk = 0.0;
		double avg_sprint_uncertainty = 0.0;
		int nn = 0;

		// Variables xij
		for (j=0; j<n; j++)
		{
			if (LP.xt[k]>0.001)
			{
				ptot += p[j]*run[j];
				objtot += LP.obj[k];
				//printf("%d ",j);
				printf("%d (%.1lf) ",j,pr[j]);
				sprint_utilization += pr[j];
				printf("%d (%.1lf) ",j,pr[j]);
				avg_sprint_risk += rcr[j];
				avg_sprint_uncertainty += run[j];
				halftotal_utility -= u[j];
				nn++;
			}
			k++;
		}
		avg_sprint_risk /= (double)(nn>0?nn:1);
		avg_sprint_uncertainty /= (double)(nn>0?nn:1);
		sprint_risks[i] = avg_sprint_risk;
		sprint_uncertainties[i] = avg_sprint_uncertainty;
		if(sprint_utilization > 0.0){
			PercentageUtilization += sprint_utilization / pmax[i];
			count++;
			NSprintsUsed++;
		}
		if(halftotal_utility <= 0.0 && HalfUtilitiSprint < 0){
			HalfUtilitiSprint = NSprintsUsed;
		}

		// Variables yij
		for (j=0; j<n; j++)
		{
			if (YMap[i*n+j]<0)
				continue;

			if (LP.xt[kk]>0.001)
			{
				objtot += LP.obj[kk];
		//		printf("y(%d,%d) = %fl   objcoef = %lf\n",i,j,LP.xt[kk],LP.obj[kk]);
			}
			kk++;
		}

		objval += objtot;
		printf("\n  ptot=%lf     objtot=%lf\n",ptot,objtot);

	}

	double mean_risk = 0.0;
	double mean_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		mean_risk += sprint_risks[i];
		mean_uncertainty += sprint_uncertainties[i];
	}

	if (NSprintsUsed > 0) {
		mean_risk /= NSprintsUsed;
		mean_uncertainty /= NSprintsUsed;
	}

	double var_risk = 0.0;
	double var_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		double dr = sprint_risks[i] - mean_risk;
		double du = sprint_uncertainties[i] - mean_uncertainty;

		var_risk += dr * dr;
		var_uncertainty += du * du;
	}

	if (NSprintsUsed > 0) {
		DeviationRisk = sqrt(var_risk / NSprintsUsed);
		DeviationUncertaintyRisk = sqrt(var_uncertainty / NSprintsUsed);
	} else {
		DeviationRisk = 0.0;
		DeviationUncertaintyRisk = 0.0;
	}

	delete [] sprint_risks;
	delete [] sprint_uncertainties;

	PercentageUtilization /= (double)count;

	printf("\n Number of Occupied Sprints: %ld", NSprintsUsed);
	printf("\n Average Sprint Utilization: %lf", PercentageUtilization);
	printf("\n Sprint Risk Standard Deviation: %lf", DeviationRisk);
	printf("\n Sprint Uncertainty Standard Deviation: %lf", DeviationUncertaintyRisk);
	if (HalfUtilitiSprint >= 0)
		printf("\n Sprint where Half of Total Utility is achieved: %ld", HalfUtilitiSprint);
	else
		printf("\n Half of Total Utility is not achieved in any sprint.");

	printf("\n Zopt=%lf\n",objval);

	// Free auxiliary data structures
	delete [] YMap;

	return 0;
}


//-----------------------------------------------------------------------------
//  Optimize heuristically the Agile Schedule                                    
//-----------------------------------------------------------------------------
//
int AgileOpt::OptimizeHeu(void)
{
	long i,j,k,h,j1,kk;
	long kount;
	double obj;
	double ptot;
	double objtot;
	double objval;
	long nYMap;
	long *YMap;
	long *Flag;

	double zheu;

	// Initialize the emerging heuristic solution
	Zheu = 0.0;
	objval = 0.0;
	printf("\n\nSoluzione: \n\n");

	// Auxiliary array to map variables "yij"
	Flag = new long [n];
	YMap = new long [n];
	nYMap = 0;
	for (j=0; j<n; j++)
	{
		Flag[j] = 0;
		if ((nY[j]>1)&&(a[j]>Prec))
		//if (a[j]>-Inf)
		{
			YMap[j]=nYMap;
			nYMap++;
		}
		else
			YMap[j]=-1;
	}


	// Initialize SprintAssign matrix (m x n)
	int** SprintAssign = new int* [m];
	for (i = 0; i < m; i++)
	{
		SprintAssign[i] = new int[n];
		for (j = 0; j < n; j++)
			SprintAssign[i][j] = 0;
	}

	CplexObj LP;

	for (i=0; i<m; i++)
	{
		// LP local into the loop body

		// Set problem size
		LP.ncols = n+nYMap;
		LP.nrows = 1 + 2*nYMap;  // OR/AND constraints are added later
		LP.nz = (n)*(2+n) + (nYMap)*(1+n);  // We overestimate |Yj|=n
		if (cfg.Sentinel)
		{
			LP.ncols += 1;
			LP.nrows += 1;
			LP.nz += n+1;
		}
		LP.MallocCols();

		// Load Columns/Variables xj in the Cplex data structure
		k = 0;
		kount = 0;
		for (j=0; j<n; j++)
		{
			//k = i*n+j;
			LP.obj[k] = (original_m-i)*u[j]*rcr[j];  // or ur[j]
			// LP.obj[k] = u[j]*rcr[j];

			LP.matbeg[k] = kount;

			// Knapsack Constraint
			LP.matind[kount] = 0;
			// LP.matval[kount] = p[j]*run[j]; // or pr[j]
			LP.matval[kount] = pr[j];
			kount++;

			// Linking Constraints xj and yj
			for (h=0; h<nY[j]; h++)
			{
				j1 = Y[j][h];
				if ((YMap[j1]>=0)&&(j1!=j))
				{
					LP.matind[kount] = 1+YMap[j1];
					LP.matval[kount] = -1.;
					kount++;
				}
			}

			if (YMap[j]>=0)
			{
				LP.matind[kount] = 1+nYMap+YMap[j];
				LP.matval[kount] = -(double)(nY[j]-1);
				kount++;
			}

			// Sentinel Constraint
			if (cfg.Sentinel)
			{
				LP.matind[kount] = 1 + 2*nYMap;
				LP.matval[kount] = 1.0;
				kount++;
			}

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;  // prima 1
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'B';

			LP.lb[k] = 0.0;
			if (Flag[j]==0)
				LP.ub[k] = 1.0;
			else
				LP.ub[k] = 0.0;
			k++;
		}

		// Load Columns/Variables yij in the Cplex data structure
		kk = 0;
		for (j=0; j<n; j++)
		{
			if (YMap[kk]<0)
			{
				kk++;
				continue;
			}

			//k = i*n+j;
			if (nY[j]<2)
				LP.obj[k] = 0.0;
			else
				LP.obj[k] = (original_m-i)*u[j]*a[j]/(nY[j]-1);
				// LP.obj[k] = a[j]/(nY[j]-1);

			LP.matbeg[k] = kount;

			LP.matind[kount] = 1+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matind[kount] = 1+nYMap+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'C';
			//LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		if (cfg.Sentinel)
		{
			// Load Columns/Variables related to Sentinel constraints in the Cplex data structure
			kk = 0;

			LP.obj[k] = 0.0;

			LP.matbeg[k] = kount;

			LP.matind[kount] = 1 + 2*nYMap;
			LP.matval[kount] = -1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 1000; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		LP.matbeg[k]=kount;

		// Load Constraints in the Cplex data structure
		kount=0;
		LP.rhs[kount]=pmax[i];
		LP.sense[kount]='L';
		kount++;

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		if (cfg.Sentinel)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='E';
			kount++;
		}

		// Maximization problem
		LP.minmax=-1;

		// Load Problem
		LP.CopyLP();

		// Set MIP
		LP.SetMIP(cfg.TimeLimit);

		// Setup Custom Cutting Plane
		LP.CuttingPlaneHeu(n,m,ur,pr,pmax,nY,Y,nUOR,UOR,nUAND,UAND,Flag,FDepP,cfg.Pred,cfg.Cover,cfg.Lifting,cfg.KnapSol,cfg.MaxCuts,&Cuts,0);

		// Optimize
		LP.SolveMIP(&zheu,&Gap,&Nodes,&Cuts);

		Zheu += zheu;

		// Read the solution
		ptot=0.0;
		objtot=0.0;
		printf("Sprint %d (pmax=%lf): ",i,pmax[i]);

		// Variables xj
		for (j=0; j<n; j++)
		{
			if (LP.xt[j]>0.001)
			{
				Flag[j] = 1;
				ptot += p[j]*run[j];
				objtot += LP.obj[j];
				//printf("%d ",j);
				printf("%d (%.1lf) ",j,pr[j]);
				SprintAssign[i][j] = 1;
			}
		}

		// Variables yij
		kk=n;
		for (j=0; j<n; j++)
		{
			if (YMap[j]<0)
				continue;

			if (LP.xt[kk]>0.001)
			{
				objtot += LP.obj[kk];
			}
			kk++;
		}

		objval += objtot;
		printf("\n  ptot=%lf    objtot=%lf\n",ptot,objtot);

	}

	
	NSprintsUsed =0;
	PercentageUtilization = 0;
	DeviationRisk = 0.0;
	DeviationUncertaintyRisk = 0.0;
	HalfUtilitiSprint = -1;
	double *sprint_risks = new double[m];
	double *sprint_uncertainties = new double[m];
	for (i = 0; i < m; i++) {
		sprint_risks[i] = 0.0;
		sprint_uncertainties[i] = 0.0;
	}

	double halftotal_utility = 0.0;
	for (j=0; j<n; j++)
		halftotal_utility += u[j];
	halftotal_utility /= 2.0;

	int count = 0;
	for (i=0; i<m; i++)
	{
		double sprint_utilization = 0.0;

		double avg_sprint_risk = 0.0;
		double avg_sprint_uncertainty = 0.0;
		int nn = 0;

		// Variables xij
		for (j=0; j<n; j++)
		{
			if (SprintAssign[i][j])
			{
				sprint_utilization += pr[j];
				printf("%d (%.1lf) ",j,pr[j]);
				avg_sprint_risk += rcr[j];
				avg_sprint_uncertainty += run[j];
				halftotal_utility -= u[j];
				nn++;
			}
			k++;
		}
		avg_sprint_risk /= (double)(nn>0?nn:1);
		avg_sprint_uncertainty /= (double)(nn>0?nn:1);
		sprint_risks[i] = avg_sprint_risk;
		sprint_uncertainties[i] = avg_sprint_uncertainty;
		if(sprint_utilization > 0.0){
			PercentageUtilization += sprint_utilization / pmax[i];
			count++;
			NSprintsUsed++;
		}
		if(halftotal_utility <= 0.0 && HalfUtilitiSprint < 0){
			HalfUtilitiSprint = NSprintsUsed;
		}

	}

	double mean_risk = 0.0;
	double mean_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		mean_risk += sprint_risks[i];
		mean_uncertainty += sprint_uncertainties[i];
	}

	if (NSprintsUsed > 0) {
		mean_risk /= NSprintsUsed;
		mean_uncertainty /= NSprintsUsed;
	}

	double var_risk = 0.0;
	double var_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		double dr = sprint_risks[i] - mean_risk;
		double du = sprint_uncertainties[i] - mean_uncertainty;

		var_risk += dr * dr;
		var_uncertainty += du * du;
	}

	if (NSprintsUsed > 0) {
		DeviationRisk = sqrt(var_risk / NSprintsUsed);
		DeviationUncertaintyRisk = sqrt(var_uncertainty / NSprintsUsed);
	} else {
		DeviationRisk = 0.0;
		DeviationUncertaintyRisk = 0.0;
	}

	delete [] sprint_risks;
	delete [] sprint_uncertainties;

	PercentageUtilization /= (double)count;

	printf("\n Number of Occupied Sprints: %ld", NSprintsUsed);
	printf("\n Average Sprint Utilization: %lf", PercentageUtilization);
	printf("\n Sprint Risk Standard Deviation: %lf", DeviationRisk);
	printf("\n Sprint Uncertainty Standard Deviation: %lf", DeviationUncertaintyRisk);
	if (HalfUtilitiSprint >= 0)
		printf("\n Sprint where Half of Total Utility is achieved: %ld", HalfUtilitiSprint);
	else
		printf("\n Half of Total Utility is not achieved in any sprint.");

	printf("\n  Zheu=%lf\n",Zheu);

	// Free auxiliary data structures
	delete [] Flag;
	delete [] YMap;

	return 0;
}

int AgileOpt::OptimizeHeu_Improved(void)
{
	long i, j, k, h, j1, kk;
	long kount;
	double obj;
	double ptot;
	double objtot;
	double objval;
	long nYMap;
	long* YMap;
	long* Flag;

	double zheu;

	// Initialize the emerging heuristic solution
	Zheu = 0.0;
	objval = 0.0;
	printf("\n\nSoluzione: \n\n");

	// Auxiliary array to map variables "yij"
	Flag = new long[n];
	YMap = new long[n];
	nYMap = 0;
	for (j = 0; j < n; j++)
	{
		Flag[j] = 0;
		if ((nY[j] > 1) && (a[j] > Prec))
		{
			YMap[j] = nYMap;
			nYMap++;
		}
		else
			YMap[j] = -1;
	}

	// Initialize SprintAssign matrix (m x n)
	int** SprintAssign = new int* [m];
	for (i = 0; i < m; i++)
	{
		SprintAssign[i] = new int[n];
		for (j = 0; j < n; j++)
			SprintAssign[i][j] = 0;
	}

	// LP local into the loop body
	CplexObj LP;

	// Set problem size
	LP.ncols = n + nYMap;
	LP.nrows = 1 + 2 * nYMap;  // OR/AND constraints are added later
	LP.nz = (n) * (2 + n) + (nYMap) * (1 + n);  // We overestimate |Yj|=n
	if (cfg.Sentinel)
	{
		LP.ncols += 1;
		LP.nrows += 1;
		LP.nz += n + 1;
	}
	LP.MallocCols();

	for (i = 0; i < m; i++)
	{
		// --- COSTRUZIONE DEL MIP ---
		k = 0;
		kount = 0;
		for (j = 0; j < n; j++)
		{
			// Funzione obiettivo migliorata
			double baseValue = u[j] * rcr[j];
			double earlyFactor = (double)(m - i);

			// Verifica dipendenze pronte
			int depsReady = 0;
			double ratio = 0.0;
			if (nY == 0) {
				ratio = 1.0;
			}
			else {
				for (int h = 0; h < nY[j]; h++) {
					int dep = Y[j][h];
					if (Flag[dep] == 1) {
						depsReady++;
					}
				}
				ratio = (double)depsReady / nY[j];
			}
			double depPenalty = 0.25 + 0.75 * ratio; 

			// Conta quante storie j sblocca
			int unlocks = 0;
			for (int k2 = 0; k2 < n; k2++) {
				for (int h = 0; h < nY[k2]; h++) {
					if (Y[k2][h] == j)
						unlocks++;
				}
			}
			double unlockBonus = 0.5 * unlocks;

			LP.obj[k] = earlyFactor * baseValue * depPenalty + unlockBonus;

			LP.matbeg[k] = kount;

			// Knapsack Constraint
			LP.matind[kount] = 0;
			LP.matval[kount] = pr[j];
			kount++;

			// Linking Constraints xj and yj
			for (h = 0; h < nY[j]; h++)
			{
				j1 = Y[j][h];
				if ((YMap[j1] >= 0) && (j1 != j))
				{
					LP.matind[kount] = 1 + YMap[j1];
					LP.matval[kount] = -1.;
					kount++;
				}
			}

			if (YMap[j] >= 0)
			{
				LP.matind[kount] = 1 + nYMap + YMap[j];
				LP.matval[kount] = -(double)(nY[j] - 1);
				kount++;
			}

			// Sentinel Constraint
			if (cfg.Sentinel)
			{
				LP.matind[kount] = 1 + 2 * nYMap;
				LP.matval[kount] = 1.0;
				kount++;
			}

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'B';

			LP.lb[k] = 0.0;
			LP.ub[k] = (Flag[j] == 0) ? 1.0 : 0.0;
			k++;
		}

		// Columns/Variables yij
		kk = 0;
		for (j = 0; j < n; j++)
		{
			if (YMap[kk] < 0)
			{
				kk++;
				continue;
			}

			if (nY[j] < 2)
				LP.obj[k] = 0.0;
			else
				LP.obj[k] = (original_m - i) * u[j] * a[j] / (nY[j] - 1);

			LP.matbeg[k] = kount;

			LP.matind[kount] = 1 + YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matind[kount] = 1 + nYMap + YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'C';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		// Sentinel columns
		if (cfg.Sentinel)
		{
			kk = 0;
			LP.obj[k] = 0.0;
			LP.matbeg[k] = kount;
			LP.matind[kount] = 1 + 2 * nYMap;
			LP.matval[kount] = -1.0;
			kount++;
			LP.matcnt[k] = kount - LP.matbeg[k];
			LP.indices[k] = k;
			LP.priority[k] = 1000;
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'I';
			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		LP.matbeg[k] = kount;

		// Constraints
		kount = 0;
		LP.rhs[kount] = pmax[i];
		LP.sense[kount] = 'L';
		kount++;

		for (k = 0; k < nYMap; k++)
		{
			LP.rhs[kount] = 0.0;
			LP.sense[kount] = 'L';
			kount++;
		}
		for (k = 0; k < nYMap; k++)
		{
			LP.rhs[kount] = 0.0;
			LP.sense[kount] = 'L';
			kount++;
		}
		if (cfg.Sentinel)
		{
			LP.rhs[kount] = 0.0;
			LP.sense[kount] = 'E';
			kount++;
		}

		LP.minmax = -1;

		LP.CopyLP();
		LP.SetMIP(cfg.TimeLimit);

		// Cutting Plane
		LP.CuttingPlaneHeu_new(n, m, ur, pr, pmax, nY, Y, nUOR, UOR, nUAND, UAND, Flag, FDepP, cfg.Pred, cfg.Cover, cfg.Lifting, cfg.KnapSol, cfg.MaxCuts, &Cuts, 0);

		LP.SolveMIP(&zheu, &Gap, &Nodes, &Cuts);

		// --- Lettura soluzione ---
		ptot = 0.0;
		objtot = 0.0;
		printf("Sprint %d (pmax=%lf): ", i, pmax[i]);

		// Registrazione assegnazioni
		for (j = 0; j < n; j++)
		{
			if (LP.xt[j] > 0.001)
			{
				Flag[j] = 1;
				SprintAssign[i][j] = 1;
			}
			else
			{
				SprintAssign[i][j] = 0;
			}
		}

		// --- Calcolo Zheu aggiornato ---
		Zheu = 0.0;
		for (long si = 0; si <= i; si++)
			for (long sj = 0; sj < n; sj++)
				if (SprintAssign[si][sj])
				{
					Zheu += (original_m - si) * u[sj] * rcr[sj];
					if (nY[sj] > 1) Zheu += a[sj] / (nY[sj] - 1);
				}

		// Output finale dello sprint
		ptot = 0.0;
		objtot = 0.0;
		for (j = 0; j < n; j++)
		{
			if (SprintAssign[i][j])
			{
				ptot += p[j] * run[j];
				objtot += LP.obj[j];
				printf("%d (%.1lf) ", j, pr[j]);
			}
		}
		printf("\n  ptot=%lf    objtot=%lf\n", ptot, objtot);
	}


	NSprintsUsed =0;
	PercentageUtilization = 0;
	DeviationRisk = 0.0;
	DeviationUncertaintyRisk = 0.0;
	HalfUtilitiSprint = -1;
	double *sprint_risks = new double[m];
	double *sprint_uncertainties = new double[m];
	for (i = 0; i < m; i++) {
		sprint_risks[i] = 0.0;
		sprint_uncertainties[i] = 0.0;
	}

	double halftotal_utility = 0.0;
	for (j=0; j<n; j++)
		halftotal_utility += u[j];
	halftotal_utility /= 2.0;

	int count = 0;
	for (i=0; i<m; i++)
	{
		double sprint_utilization = 0.0;

		double avg_sprint_risk = 0.0;
		double avg_sprint_uncertainty = 0.0;
		int nn = 0;

		// Variables xij
		for (j=0; j<n; j++)
		{
			if (SprintAssign[i][j])
			{
				sprint_utilization += pr[j];
				printf("%d (%.1lf) ",j,pr[j]);
				avg_sprint_risk += rcr[j];
				avg_sprint_uncertainty += run[j];
				halftotal_utility -= u[j];
				nn++;
			}
			k++;
		}
		printf("\n");
		printf("Sprint Utilization: %lf\n", sprint_utilization);
		printf("\n");
		avg_sprint_risk /= (double)(nn>0?nn:1);
		avg_sprint_uncertainty /= (double)(nn>0?nn:1);
		sprint_risks[i] = avg_sprint_risk;
		sprint_uncertainties[i] = avg_sprint_uncertainty;
		if(sprint_utilization > 0.0){
			PercentageUtilization += sprint_utilization / pmax[i];
			count++;
			NSprintsUsed++;
		}
		if(halftotal_utility <= 0.0 && HalfUtilitiSprint < 0){
			HalfUtilitiSprint = NSprintsUsed;
		}

	}

	double mean_risk = 0.0;
	double mean_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		mean_risk += sprint_risks[i];
		mean_uncertainty += sprint_uncertainties[i];
	}

	if (NSprintsUsed > 0) {
		mean_risk /= NSprintsUsed;
		mean_uncertainty /= NSprintsUsed;
	}

	double var_risk = 0.0;
	double var_uncertainty = 0.0;

	for (i = 0; i < NSprintsUsed; i++) {
		double dr = sprint_risks[i] - mean_risk;
		double du = sprint_uncertainties[i] - mean_uncertainty;

		var_risk += dr * dr;
		var_uncertainty += du * du;
	}

	if (NSprintsUsed > 0) {
		DeviationRisk = sqrt(var_risk / NSprintsUsed);
		DeviationUncertaintyRisk = sqrt(var_uncertainty / NSprintsUsed);
	} else {
		DeviationRisk = 0.0;
		DeviationUncertaintyRisk = 0.0;
	}

	delete [] sprint_risks;
	delete [] sprint_uncertainties;

	PercentageUtilization /= (double)count;

	printf("\n Number of Occupied Sprints: %ld", NSprintsUsed);
	printf("\n Average Sprint Utilization: %lf", PercentageUtilization);
	printf("\n Sprint Risk Standard Deviation: %lf", DeviationRisk);
	printf("\n Sprint Uncertainty Standard Deviation: %lf", DeviationUncertaintyRisk);
	if (HalfUtilitiSprint >= 0)
		printf("\n Sprint where Half of Total Utility is achieved: %ld", HalfUtilitiSprint);
	else
		printf("\n Half of Total Utility is not achieved in any sprint.");

	printf("\n  Zheu=%lf\n", Zheu);

	// Free auxiliary data structures
	delete[] Flag;
	delete[] YMap;
	for (i = 0; i < m; i++) delete[] SprintAssign[i];
	delete[] SprintAssign;

	return 0;
}


//-----------------------------------------------------------------------------
//  Optimize heuristically the Agile Schedule by a Lagrangian Heuristc                                   
//-----------------------------------------------------------------------------
//
int AgileOpt::OptimizeLagrHeu(void)
{
	long i,j,k,h,j1,kk;
	long it;
	long kount;
	double obj;
	double ptot;
	long nYMap;
	long *YMap;
	long *Flag;
	double *urpen;
	double *apen;
	double *xlr;
	double *ylr;
	double *xheu;
	double *yheu;
	double *Lambda;
	double *LambdaOR;
	double *LambdaAND;
	double *LambdaY1;
	double *LambdaY2;
	double *subgr;
	double *subgrOR;
	double *subgrAND;
	double *subgrY1;
	double *subgrY2;

	double Konst;
	double zheut;
	double zheu;
	double zlr;
	double bestzheu;
	double bestzlr;
	double alpha;
	double zgap;
	double mingap;
	long maxiter;
	long maxitgap;
	long maxbadit;
	long badit;

	int status,bol,bolheu;
	long t1,t2;
	double dt;

	t1 = clock();

//2012	alpha = 0.1;  // 0.1 (0.2 converge piu' veloce, pero' 0.1 e' piu' stabile)
	alpha = 30.0;  // Best=24-30   (24) (10) 0.1 (0.2 converge piu' veloce, pero' 0.1 e' piu' stabile)
	maxiter = 5000;
	maxitgap = 50; // Best=50
	maxbadit = 10;  // Best=12-10  (12) (15) 10 (20 converge piu' veloce)
	mingap = 0.01; // Best=0.01
	bestzheu = 0.0;
	bestzlr = +1000000000.;
	zgap = bestzlr;

	// Initialize the emerging heuristic solution
	Zheu = 0.0;
	zheu = 0.0;
	Zlagr = +1000000000.;
	printf("\n\nSoluzione: \n\n");

	// Try with the greedy procedure (one shot)
	OptimizeHeu();
	zheu = Zheu;
	bestzheu = Zheu;

	// Allocate Lagrangian Penalties
	urpen = new double [n*m];
	apen = new double [n*m];
	xlr = new double [n*m];
	ylr = new double [n*m];
	xheu = new double [n*m];
	yheu = new double [n*m];
	Lambda = new double [n];
	LambdaOR = new double [n*m];
	LambdaAND = new double [n*m];
	LambdaY1 = new double [n*m];
	LambdaY2 = new double [n*m];
	subgr = new double [n];
	subgrOR = new double [n*m];
	subgrAND = new double [n*m];
	subgrY1 = new double [n*m];
	subgrY2 = new double [n*m];
	for (j=0; j<n; j++)
	{
		Lambda[j] = 0.0;
		subgr[j] = 0.0;
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] = (m-i)*ur[j];
			if (nY[j]<2)
				apen[k] = 0.0;
			else
				apen[k] = (m-i)*u[j]*a[j]/(nY[j]-1);
			LambdaOR[k] = 0.0;
			LambdaAND[k] = 0.0;
			LambdaY1[k] = 0.0;
			LambdaY2[k] = 0.0;
			subgrOR[k] = 0.0;
			subgrAND[k] = 0.0;
			subgrY1[k] = 0.0;
			subgrY2[k] = 0.0;
			xlr[k] = 0.0;
			ylr[k] = 0.0;
			xheu[k] = 0.0;
		}
	}

	// Auxiliary array to map variables "yij"
	Flag = new long [n];
	YMap = new long [n];
	nYMap = 0;
	for (j=0; j<n; j++)
	{
		Flag[j] = 0;
		if ((nY[j]>1)&&(a[j]>Prec))
		//if (a[j]>-Inf)
		{
			YMap[j]=nYMap;
			nYMap++;
		}
		else
			YMap[j]=-1;
	}

	for (it=0; it<1000; it++)
	{
		// Compute the penalized price 
		//ComputePenPrice(urpen,Lambda,LambdaOR,LambdaAND,&Konst);
		ComputePenPriceDP(urpen,apen,Lambda,LambdaOR,LambdaAND,LambdaY1,LambdaY2,&Konst);

		// Compute the Lagrangian Problem Solution
		//SolveLR(urpen,apen,nYMap,YMap,Flag,Konst,&zlr,xlr,ylr);
		SolveLRDP(urpen,apen,nYMap,YMap,Flag,Konst,&zlr,xlr,ylr);
		if (bestzlr>zlr) 
		{
			bolheu = 1;
			bestzlr = zlr;
			badit = 0;
		}
		else
		{
			bolheu = 0;
			badit++;
			if (badit>=maxbadit)
			{
				alpha = 0.85*alpha;
				badit = 0;
			}
		}

		// Compute an feasible solution using a Lagrangian Heuristic
		//LagrHeu(urpen,apen,nYMap,YMap,Flag,&zheu,xheu,yheu); 
		//if (bestzheu<zheu) 
		//	bestzheu = zheu;
		if ((bolheu)||(it<0))
		{
			//LagrHeuDP(urpen,nYMap,YMap,Flag,&zheu,xheu,yheu); 
			LagrHeuDPBack(urpen,nYMap,YMap,Flag,&zheu,xheu,yheu); 
			if (zheu>0.0)
			{
				// HeuSwaps(&zheu,xheu,yheu);
				if (bestzheu<zheu) 
					bestzheu = zheu;
				if (cfg.HeuBest*bestzheu<zheu) 
				{
					//if (bestzheu<zheu) 
					//	bestzheu = zheu;
					LagrHeu(urpen,apen,nYMap,YMap,Flag,&zheu,xheu,yheu); 
					if (zheu>0.0)
					{
						//HeuSwaps(&zheu,xheu,yheu);
						if (bestzheu<zheu) 
							bestzheu = zheu;
					}
				}
			}
		}

		//LagrHeuDP_Old(urpen,nYMap,YMap,Flag,&zheu,xheu); 
		//if (bestzheu<zheu) 
		//	bestzheu = zheu;

		if ((it%maxitgap)==0)
		{
			if ((100.*(zgap-bestzlr)/bestzlr)<mingap)
				break;
			else
				zgap = bestzlr;
		}

		// Update penalties
		//UpdatePen(alpha,bestzheu,zlr,xlr,Lambda,LambdaOR,LambdaAND,subgr,subgrOR,subgrAND);
		UpdatePenDP(alpha,bestzheu,zlr,xlr,ylr,Lambda,LambdaOR,LambdaAND,LambdaY1,LambdaY2,subgr,subgrOR,subgrAND,subgrY1,subgrY2);

		if ((it%100)==0)
			printf("%d) Zlr=%.1lf    Zheu=%.1lf    BestZlr=%.1lf    BestZheu=%.1lf    Alpha=%.5lf    Konst=%.5lf\n",it,zlr,zheu,bestzlr,bestzheu,alpha,Konst);

		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		// Check Time Limit
		if (dt>cfg.TimeLimit)
			break;
	}

	Zheu = bestzheu;
	Zlagr = bestzlr;

	printf("\n  Zheu=%lf\n",Zheu);

	// Free Penalty vectors
	delete [] urpen;
	delete [] apen;
	delete [] xlr;
	delete [] ylr;
	delete [] xheu;
	delete [] yheu;
	delete [] Lambda;
	delete [] LambdaOR;
	delete [] LambdaAND;
	delete [] LambdaY1;
	delete [] LambdaY2;
	delete [] subgr;
	delete [] subgrOR;
	delete [] subgrAND;
	delete [] subgrY1;
	delete [] subgrY2;

	// Free auxiliary data structures
	delete [] Flag;
	delete [] YMap;

	return 0;
}


//-----------------------------------------------------------------------------
//  Optimize heuristically the Agile Schedule by a Lagrangian Heuristc                                   
//-----------------------------------------------------------------------------
//
int AgileOpt::ComputePenPrice(double *urpen, double *Lambda, double *LambdaOR, double *LambdaAND, double *Konst)
{
	long i,j,k;
	long i1,j1,k1;
	long jj;
	double PKonst;

	PKonst = 0.0;
	
	// Assignment constraints
	for (j=0; j<n; j++)
	{
		PKonst += Lambda[j];
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] = (m-i)*ur[j] - Lambda[j];
		}
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] += LambdaOR[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUOR[j]; jj++)
				{
					j1 = UOR[j][jj];
					k1 = i1*n+j1;
					urpen[k1] -= LambdaOR[k];
				}
			}
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] += (double)nUAND[j]*LambdaAND[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUAND[j]; jj++)
				{
					j1 = UAND[j][jj];
					k1 = i1*n+j1;
					urpen[k1] -= LambdaAND[k];
				}
			}
		}
	}

	*Konst = PKonst;

	return 0;
}


//-----------------------------------------------------------------------------
//  Optimize heuristically the Agile Schedule by a Lagrangian Heuristc                                   
//-----------------------------------------------------------------------------
//
int AgileOpt::ComputePenPriceDP(double *urpen, double *apen, double *Lambda, double *LambdaOR, double *LambdaAND,  
	                            double *LambdaY1, double *LambdaY2, double *Konst)
{
	long i,j,k;
	long i1,j1,k1;
	long jj;
	double PKonst;

	PKonst = 0.0;
	
	// Assignment constraints
	for (j=0; j<n; j++)
	{
		PKonst += Lambda[j];
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] = (m-i)*ur[j] - Lambda[j];
		}
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] += LambdaOR[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUOR[j]; jj++)
				{
					j1 = UOR[j][jj];
					k1 = i1*n+j1;
					urpen[k1] -= LambdaOR[k];
				}
			}
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			urpen[k] += (double)nUAND[j]*LambdaAND[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUAND[j]; jj++)
				{
					j1 = UAND[j][jj];
					k1 = i1*n+j1;
					urpen[k1] -= LambdaAND[k];
				}
			}
		}
	}

	// Constraints Y1
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			apen[k] = (m-i)*u[j]*a[j]/(nY[j]-1) + LambdaY1[k];
			for (jj=0; jj<nY[j]; jj++)
			{
				j1 = Y[j][jj];
				if (j1==j)
					continue;

				k1 = i*n+j1;
				urpen[k1] -= LambdaY1[k];
			}
		}
	}

	// Constraints Y2
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			apen[k] += LambdaY2[k];
			urpen[k] -= (double)(nY[j]-1)*LambdaY2[k];
		}
	}

	*Konst = PKonst;

	return 0;
}


//-----------------------------------------------------------------------------
// Compute the Lagrangian Problem Solution
//-----------------------------------------------------------------------------
//
int AgileOpt::SolveLR(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double Konst, double *zlr, double *xlr, double *ylr)
{
	long i,j,k,h,j1,kk,k1;
	long kount;
	double ptot;
	double objtot;
	double objval;
	double zlrt;

	for (j=0; j<n; j++)
	{
		Flag[j] = 0;
	}

	// Initialize the solution value with the constant contribution of the Lagrangian penalties
	objval = 0.0;
	(*zlr) = Konst;

	for (i=0; i<m; i++)
	{
		// LP local into the loop body
		CplexObj LP;

		// Set problem size
		LP.ncols = n+nYMap;
		LP.nrows = 1 + 2*nYMap;  // OR/AND constraints are added later
		LP.nz = (n)*(2+n) + (nYMap)*(1+n);  // We overestimate |Yj|=n
		if (cfg.Sentinel)
		{
			LP.ncols += 1;
			LP.nrows += 1;
			LP.nz += n+1;
		}
		LP.MallocCols();

		// Load Columns/Variables xj in the Cplex data structure
		k = 0;
		kount = 0;
		for (j=0; j<n; j++)
		{
			kk = i*n+j;
			//LP.obj[k] = (m-i)*urpen[k];
			LP.obj[k] = urpen[kk];

			LP.matbeg[k] = kount;

			// Knapsack Constraint
			LP.matind[kount] = 0;
			// LP.matval[kount] = p[j]*run[j]; // or pr[j]
			LP.matval[kount] = pr[j];
			kount++;

			// Linking Constraints xj and yj
			for (h=0; h<nY[j]; h++)
			{
				j1 = Y[j][h];
				if ((YMap[j1]>=0)&&(j1!=j))
				{
					LP.matind[kount] = 1+YMap[j1];
					LP.matval[kount] = -1.;
					kount++;
				}
			}

			if (YMap[j]>=0)
			{
				LP.matind[kount] = 1+nYMap+YMap[j];
				LP.matval[kount] = -(double)(nY[j]-1);
				kount++;
			}

			// Sentinel Constraint
			if (cfg.Sentinel)
			{
				LP.matind[kount] = 1 + 2*nYMap;
				LP.matval[kount] = 1.0;
				kount++;
			}

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;  // prima 1
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'B';

			LP.lb[k] = 0.0;
			if (Flag[j]==0)
				LP.ub[k] = 1.0;
			else
				LP.ub[k] = 0.0;
			k++;
		}

		// Load Columns/Variables yij in the Cplex data structure
		kk = 0;
		for (j=0; j<n; j++)
		{
			if (YMap[kk]<0)
			{
				kk++;
				continue;
			}

			k1 = i*n+j;
			if (nY[j]<2)
				LP.obj[k] = 0.0;
			else
				LP.obj[k] = apen[k1];
				//LP.obj[k] = (m-i)*u[j]*a[j]/(nY[j]-1);
				// LP.obj[k] = a[j]/(nY[j]-1);

			LP.matbeg[k] = kount;

			LP.matind[kount] = 1+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matind[kount] = 1+nYMap+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'C';
			//LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		if (cfg.Sentinel)
		{
			// Load Columns/Variables related to Sentinel constraints in the Cplex data structure
			kk = 0;

			LP.obj[k] = 0.0;

			LP.matbeg[k] = kount;

			LP.matind[kount] = m + n + 2*nYMap;
			LP.matval[kount] = -1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 1000; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		LP.matbeg[k]=kount;

		// Load Constraints in the Cplex data structure
		kount=0;
		LP.rhs[kount]=pmax[i];
		LP.sense[kount]='L';
		kount++;

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		if (cfg.Sentinel)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='E';
			kount++;
		}

		// Maximization problem
		LP.minmax=-1;

		// Load Problem
		LP.CopyLP();

		// Set MIP
		LP.SetMIP(cfg.TimeLimit);

		// Setup Custom Cutting Plane
		LP.CuttingPlaneHeu(n,m,ur,pr,pmax,nY,Y,nUOR,UOR,nUAND,UAND,Flag,FDepP,cfg.Pred,cfg.Cover,cfg.Lifting,cfg.KnapSol,cfg.MaxCuts,&Cuts,i);
		//LP.CuttingPlaneHeu(n,m,ur,pr,pmax,nY,Y,nUOR,UOR,nUAND,UAND,Flag,FDepP,cfg.Pred,cfg.Cover,cfg.Lifting,cfg.KnapSol,cfg.MaxCuts,&Cuts,1);

		// Optimize
		LP.SolveMIP(&zlrt,&Gap,&Nodes,&Cuts);

		(*zlr) += zlrt;

		// Read the solution
		ptot=0.0;
		objtot=0.0;
		if (iprinx>0) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);

		// Variables xj
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (LP.xt[j]>0.001)
			{
				//Flag[j] = 1;
				xlr[k] = 1.0;
				ptot += p[j]*run[j];
				objtot += LP.obj[j];
				if (iprinx>0) printf("%d ",j);
			}
			else
			{
				xlr[k] = 0.0;
			}
		}

		// Variables yj
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (YMap[j]<0)
			{
				ylr[k] = 0.0;
				continue;
			}

			kk = n + YMap[j];
			if (LP.xt[kk]>0.001)
			{
				ylr[k] = LP.xt[kk];
				objtot += LP.obj[kk];
				//(*zlr) += (m-i)*u[j]*a[j]*LP.xt[kk];
			}
			else
				ylr[k] = 0.0;
		}

		objval += objtot;
		if (iprinx>0) printf("\n  ptot=%lf    objtot=%lf\n",ptot,objtot);
	}

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zlr));

	return 0;
}


//-----------------------------------------------------------------------------
// Compute the Lagrangian Problem Solution using Dynamic Programming
//-----------------------------------------------------------------------------
//
int AgileOpt::SolveLRDP(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double Konst, double *zlr, double *xlr, double *ylr)
{
	long i,j,k,h,j1,kk;
	long kount;
	double ptot;
	double objtot;
	double objval;
	double zlrt;
	long maxw;
	Knapsack KP;

	for (j=0; j<n; j++)
	{
		Flag[j] = 0;
	}

	maxw = 0;
	for (i=0; i<m; i++)
	{
		if (maxw<(int(pmax[i]*10+0.001)))
			maxw = int(pmax[i]*10+0.001);
	}

	// Allocate Knapsack data structures
	KP.Malloc(n,maxw);

	// Initialize the solution value with the constant contribution of the Lagrangian penalties
	objval = 0.0;
	(*zlr) = Konst;

	for (i=0; i<m; i++)
	{
		// LP local into the loop body
		//CplexObj LP;

		// Variables xij
		KP.n = 0;
		KP.W = int(pmax[i]*10+0.001);
		for (j=0; j<n; j++) 
		{
			k = n*i+j;
			if (urpen[k]>Prec)
			{
				KP.w[KP.n] = int(pr[j]*10+0.001);
				KP.p[KP.n] = urpen[k];
				KP.ind[KP.n] = j;
				KP.n++;
			}
		}

		KP.SolveReMap();
		zlrt = KP.zopt;

		// Read the solution
		ptot=0.0;
		objtot=0.0;
		if (iprinx>0) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			xlr[k] = 0.0;
		}
		for (j1=0;j1<KP.n;j1++)
		{
			j = KP.ind[j1];
			k = i*n+j;
			xlr[k] = KP.x[j1];
			if (xlr[k]>Prec)
			{
				ptot += p[j]*run[j];
				objtot += KP.p[j1];
				if (iprinx>0) printf("%d ",j);
			}
		}

		// Variables yij
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			ylr[k] = 0.0;
		}
		for (j=0; j<n; j++)
		{
			if (nY[j]<2)
				continue;

			k = n*i+j;
			if (apen[k]>Prec)
			{
				ylr[k] = (double)(nY[j]-1);
				zlrt += ylr[k] * apen[k];
				objtot += ylr[k] * apen[k];
			}
			else
			{
				ylr[k] = 0.0;
			}
		}

		objval += objtot; 
		if (iprinx>0) printf("\n  ptot=%lf    objtot=%lf\n",ptot,objtot);

		(*zlr) += zlrt;

	}

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zlr));

	return 0;
}


//-----------------------------------------------------------------------------
// Try to move j1 from i1 to i
//-----------------------------------------------------------------------------
//
double AgileOpt::TryMove(double *xheu, double *yheu, long j1, long i1, long i)
{
	long jj,j2,k1,k2;
	long konta;
	double gain;

	k1 = i1*n+j1;

	// Gain for moving j1 from i1 to i
	gain = (double)(i1-i)*ur[j1]; // Saving due to x[k] 
	if (nY[j1]>1)
	{
		// Saving due to y[k] 
		gain -= (double)(m-i1)*yheu[k1]*u[j1]*a[j1]/(nY[j1]-1);
		konta = 0;
		for (jj=0; jj<nY[j1]; jj++)
		{
			j2 = Y[j1][jj];
			if (j2==j1) 
				continue;
			k2 = i*n+j2;
			if (xheu[k2]>0.999)
			{
				konta++;
				gain += (double)(m-i)*(+1.)*u[j2]*a[j2]/(nY[j2]-1);
			}
			k2 = i1*n+j2;
			if (xheu[k2]>0.999)
				gain += (double)(m-i1)*(-1.)*u[j2]*a[j2]/(nY[j2]-1); 
		}
		gain += (double)((m-i)*konta)*u[j1]*a[j1]/(nY[j1]-1);
	}

	return gain;
}


//-----------------------------------------------------------------------------
// Try to move j1 from i1 to i
//-----------------------------------------------------------------------------
//
double AgileOpt::Move(double *xheu, double *yheu, long j1, long i1, long i)
{
	long jj,j2,k1,k2;
	long konta;

	k1 = i1*n+j1;

	// Moving j1 from i1 to i
	xheu[k1] = 0.0;
	yheu[k1] = 0.0;
	konta = 0;
	if (nY[j1]>1)
	{
		for (jj=0; jj<nY[j1]; jj++)
		{
			j2 = Y[j1][jj];
			if (j2==j1) 
				continue;
			k2 = i*n+j2;
			if (xheu[k2]>0.999)
			{
				konta++;
				yheu[k2] += 1.;
			}
			k2 = i1*n+j2;
			if (xheu[k2]>0.999)
				yheu[k2] -= 1.; 
		}
	}
	k2 = i*n+j1; 
	xheu[k2] = 1.0;
	yheu[k2] = (double)konta;

	return 0.0;
}


//-----------------------------------------------------------------------------
// Improve a heuristic solutions by swaping user stories
//-----------------------------------------------------------------------------
//
int AgileOpt::HeuSwaps(double *zheu, double *xheu, double *yheu)
{
	long i,j,k;
	long i1,j1,k1;
	long jj,j2,k2;
	long konta;
	long konta1;
	double zold;
	double ptot,ptot1;
	double gain;
	double zheuf;
	int status,bolf;
	int bol;

	zold = *zheu;

TryAgain:

	bol = 0;

	// Swaps 1-0
	for (i=0; i<m-1; i++)
	{
		ptot = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (xheu[k]>0.999)
				ptot += pr[j];
		}
		for (i1=i+1; i1<m; i1++)
		{
			for (j1=0; j1<n; j1++)
			{
				k1 = i1*n+j1;
				if (xheu[k1]>0.999)
				{
					if ((ptot+pr[j1])<pmax[i]+Prec)  // Can j1 be anticipated?
					{
						// Check if the objective fuction increase!
						//gain = TryMove(xheu,yheu,j1,i1,i);

						//if (gain>Prec)
						//{
						ptot += pr[j1];
						Move(xheu,yheu,j1,i1,i);

						// Check feasibility
						status = Feasibility(xheu,yheu,&zheuf,&bolf);
						gain = zheuf-(*zheu); 
						if ((status)||(gain<Prec))
						{
							// Backtracking
							ptot += -pr[j1];
							Move(xheu,yheu,j1,i,i1);
						}
						else
						{
							//(*zheu) += gain;
							(*zheu) = zheuf;
							bol = 1;
							goto Next0;
						}
						//}
					}
				}

Next0:;

			}
		}
	}

	// Swaps 1-1
	for (i=0; i<m-1; i++)
	{
		ptot = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (xheu[k]>0.999)
				ptot += pr[j];
		}

		for (i1=i+1; i1<m; i1++)
		{
			ptot1 = 0.0;
			for (j1=0; j1<n; j1++)
			{
				k1 = i1*n+j1;
				if (xheu[k1]>0.999)
					ptot1 += pr[j1];
			}

			for (j=0; j<n; j++)
			{
				k = i*n+j;
				if (xheu[k]<0.999)
					continue;

				for (j1=0; j1<n; j1++)
				{
					if (j==j1)
						continue;

					k1 = i1*n+j1;
					if (xheu[k1]>0.999)
					{ 
						// Can j and j1 be swapped?
						if (((ptot-pr[j]+pr[j1])<pmax[i]+Prec)&&((ptot1+pr[j]-pr[j1])<pmax[i]+Prec)) 
						{
							// Check if the objective fuction increase!

							// Gain for moving j1 from i1 to i
							//gain = TryMove(xheu,yheu,j1,i1,i);

							// Gain for moving j from i to i1
							//gain += TryMove(xheu,yheu,j,i,i1);

							//gain=1.;
							//if (gain>Prec)
							//{
							ptot += pr[j1]-pr[j];
							ptot1 += pr[j]-pr[j1];

							Move(xheu,yheu,j1,i1,i);
							Move(xheu,yheu,j,i,i1);

							// Check feasibility
							status = Feasibility(xheu,yheu,&zheuf,&bolf);
							gain = zheuf-(*zheu); 
							if ((status)||(gain<Prec))
							{
								// Backtracking
								ptot += pr[j]-pr[j1];
								ptot1 += pr[j1]-pr[j];

								Move(xheu,yheu,j1,i,i1);
								Move(xheu,yheu,j,i1,i);
							}
							else
							{
								//(*zheu) += gain;
								(*zheu) = zheuf;
								bol = 1;
								goto Next1;
							}
							//}
						}
					}
				}
Next1:;

			}
		}
	}

	//goto END;

	// Swaps 2-1
	for (i=0; i<m; i++)
	{
		ptot = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (xheu[k]>0.999)
				ptot += pr[j];
		}

		for (i1=0; i1<m; i1++)
		{
			if (i==i1)
				continue;

			ptot1 = 0.0;
			for (j1=0; j1<n; j1++)
			{
				k1 = i1*n+j1;
				if (xheu[k1]>0.999)
					ptot1 += pr[j1];
			}

			for (j=0; j<n-1; j++)
			{
				k = i*n+j;
				if (xheu[k]<0.999)
					continue;

				for (j1=j+1; j1<n; j1++)
				{
					k1 = i*n+j1;
					if (xheu[k1]<0.999)
						continue;

					for (j2=0; j2<n; j2++)
					{
						k2 = i1*n+j2;
						if (xheu[k2]>0.999)
						{ 
							// Can (j and j1) be swapped with j2?
							if (((ptot-pr[j]-pr[j1]+pr[j2])<pmax[i]+Prec)&&((ptot1+pr[j]+pr[j1]-pr[j2])<pmax[i]+Prec)) 
							{
								// Check if the objective fuction increase!

								//// Gain for moving j2 from i1 to i
								//gain = TryMove(xheu,yheu,j2,i1,i);

								//// Gain for moving j1 from i to i1
								//gain = TryMove(xheu,yheu,j1,i,i1);

								//// Gain for moving j from i to i1
								//gain += TryMove(xheu,yheu,j,i,i1);

								//if (gain>Prec)
								//{
								ptot += pr[j2]-pr[j]-pr[j1];
								ptot1 += pr[j]+pr[j1]-pr[j2];

								Move(xheu,yheu,j2,i1,i);
								Move(xheu,yheu,j1,i,i1);
								Move(xheu,yheu,j,i,i1);

								// Check feasibility
								status = Feasibility(xheu,yheu,&zheuf,&bolf);
								gain = zheuf-(*zheu); 
								if ((status)||(gain<Prec))
								{
									// Backtracking
									ptot += pr[j]+pr[j1]-pr[j2];
									ptot1 += pr[j2]-pr[j]-pr[j1];

									Move(xheu,yheu,j2,i,i1);
									Move(xheu,yheu,j1,i1,i);
									Move(xheu,yheu,j,i1,i);
								}
								else
								{
									//(*zheu) += gain;
									(*zheu) = zheuf;
									bol = 1;
									goto Next2;
								}
								//}
							}
						}
					}
				}

Next2:;

			}
		}
	}

	if (bol)
		goto TryAgain;

END:

	// Check feasibility
	status = Feasibility(xheu,yheu,&zheuf,&bolf);
	if (status)
	{
		//printf("ERROR...  HeuSwaps solution not feasible: status=%d\n",status);
		//fprintf(ferr,"ERROR...  HeuSwaps solution not feasible: status=%d\n",status);
		*zheu=-Inf;
		// exit(status);
	}
	else if (((*zheu)>zheuf+LPrec)||((*zheu)<zheuf-LPrec))
	{
		//printf("ERROR...  HeuSwaps value: zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		//fprintf(ferr,"ERROR...  HeuSwaps value: zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		*zheu=-Inf;
		// exit(status);
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Improve a heuristic solutions by swaping user stories
//-----------------------------------------------------------------------------
//
int AgileOpt::HeuSwapsOld(double *zheu, double *xheu, double *yheu)
{
	long i,j,k;
	long i1,j1,k1;
	long jj,j2,k2;
	long konta;
	long konta1;
	double zold;
	double ptot,ptot1;
	double gain;
	double zheuf;
	int status,bolf;

	zold = *zheu;

	// Swaps 1-0
	for (i=0; i<m-1; i++)
	{
		ptot = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (xheu[k]>0.999)
				ptot += pr[j];
		}
		for (i1=i+1; i1<m; i1++)
		{
			for (j1=0; j1<n; j1++)
			{
				k1 = i1*n+j1;
				if (xheu[k1]>0.999)
				{
					if ((ptot+pr[j1])<pmax[i]+Prec)  // Can j1 be anticipated?
					{
						// Check if the objective fuction increase!
						gain = (double)(i1-i)*ur[j1]; // Saving due to x[k] 
						if (nY[j1]>1)
						{
							// Saving due to y[k] 
							gain -= (double)(m-i1)*yheu[k1]*u[j1]*a[j1]/(nY[j1]-1);
							konta = 0;
							for (jj=0; jj<nY[j1]; jj++)
							{
								j = Y[j1][jj];
								if (j==j1) 
									continue;
								k = i*n+j;
								if (xheu[k]>0.999)
								{
									konta++;
									gain += (double)(m-i)*(+1.)*u[j]*a[j]/(nY[j]-1);
								}
								k = i1*n+j;
								if (xheu[k]>0.999)
									gain += (double)(m-i1)*(-1.)*u[j]*a[j]/(nY[j]-1); 
							}
							gain += (double)((m-i)*konta)*u[j1]*a[j1]/(nY[j1]-1);
						}
						if (gain>Prec)
						{
							(*zheu) += gain;
							ptot += pr[j1];
							xheu[k1] = 0.0;
							yheu[k1] = 0.0;
							k = i*n+j1; 
							xheu[k] = 1.0;
							yheu[k] = (double)konta;
							if (nY[j1]>1)
							{
								for (jj=0; jj<nY[j1]; jj++)
								{
									j = Y[j1][jj];
									if (j==j1) 
										continue;
									k = i*n+j;
									if (xheu[k]>0.999)
										yheu[k] += 1.;
									k = i1*n+j;
									if (xheu[k]>0.999)
										yheu[k] -= 1.; 
								}
							}
						}
					}
				}
			}
		}
	}

	// Swaps 1-1
	for (i=0; i<m-1; i++)
	{
		ptot = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (xheu[k]>0.999)
				ptot += pr[j];
		}

		for (i1=i+1; i1<m; i1++)
		{
			ptot1 = 0.0;
			for (j1=0; j1<n; j1++)
			{
				k1 = i1*n+j1;
				if (xheu[k1]>0.999)
					ptot1 += pr[j1];
			}

			for (j=0; j<n-1; j++)
			{
				k = i*n+j;
				if (xheu[k]<0.999)
					continue;

				for (j1=j+1; j1<n; j1++)
				{
					k1 = i1*n+j1;
					if (xheu[k1]>0.999)
					{ 
						// Can j and j1 be swapped?
						if (((ptot-pr[j]+pr[j1])<pmax[i]+Prec)&&((ptot1+pr[j]-pr[j1])<pmax[i]+Prec)) 
						{
							// Check if the objective fuction increase!

							// Gain for moving j1 from i1 to i
							gain = (double)(i1-i)*ur[j1]; // Saving due to x[k] 
							if (nY[j1]>1)
							{
								// Saving due to y[k] 
								gain -= (double)(m-i1)*yheu[k1]*u[j1]*a[j1]/(nY[j1]-1);
								konta1 = 0;
								for (jj=0; jj<nY[j1]; jj++)
								{
									j2 = Y[j1][jj];
									if (j2==j1) 
										continue;
									k2 = i*n+j2;
									if (xheu[k2]>0.999)
									{
										konta1++;
										gain += (double)(m-i)*(+1.)*u[j2]*a[j2]/(nY[j2]-1);
									}
									k2 = i1*n+j2;
									if (xheu[k2]>0.999)
										gain += (double)(m-i1)*(-1.)*u[j2]*a[j2]/(nY[j2]-1); 
								}
								gain += (double)((m-i)*konta1)*u[j1]*a[j1]/(nY[j1]-1);
							}

							// Gain for moving j from i to i1
							gain -= (double)(i1-i)*ur[j]; // Saving due to x[k] 
							if (nY[j]>1)
							{
								// Saving due to y[k] 
								gain -= (double)(m-i)*yheu[k]*u[j]*a[j]/(nY[j]-1);
								konta = 0;
								for (jj=0; jj<nY[j]; jj++)
								{
									j2 = Y[j][jj];
									if (j2==j) 
										continue;
									k2 = i1*n+j2;
									if (xheu[k2]>0.999)
									{
										konta++;
										gain += (double)(m-i1)*(+1.)*u[j2]*a[j2]/(nY[j2]-1);
									}
									k2 = i*n+j2;
									if (xheu[k2]>0.999)
										gain += (double)(m-i)*(-1.)*u[j2]*a[j2]/(nY[j2]-1); 
								}
								gain += (double)((m-i1)*konta)*u[j]*a[j]/(nY[j]-1);
							}
							if (gain>Prec)
							{
								(*zheu) += gain;
								ptot += pr[j1]-pr[j];
								ptot1 += pr[j]-pr[j1];

								// Moving j1 from i1 to i
								xheu[k1] = 0.0;
								yheu[k1] = 0.0;
								k2 = i*n+j1; 
								xheu[k2] = 1.0;
								yheu[k2] = (double)konta1;
								if (nY[j1]>1)
								{
									for (jj=0; jj<nY[j1]; jj++)
									{
										j2 = Y[j1][jj];
										if (j2==j1) 
											continue;
										k2 = i*n+j2;
										if (xheu[k2]>0.999)
											yheu[k2] += 1.;
										k2 = i1*n+j2;
										if (xheu[k2]>0.999)
											yheu[k2] -= 1.; 
									}
								}

								// Moving j from i to i1
								xheu[k] = 0.0;
								yheu[k] = 0.0;
								k2 = i1*n+j; 
								xheu[k2] = 1.0;
								yheu[k2] = (double)konta;
								if (nY[j]>1)
								{
									for (jj=0; jj<nY[j]; jj++)
									{
										j2 = Y[j][jj];
										if (j2==j) 
											continue;
										k2 = i1*n+j2;
										if (xheu[k2]>0.999)
											yheu[k2] += 1.;
										k2 = i*n+j2;
										if (xheu[k2]>0.999)
											yheu[k2] -= 1.; 
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// Check feasibility
	status = Feasibility(xheu,yheu,&zheuf,&bolf);
	if (status)
	{
		printf("ERROR...  HeuSwaps solution not feasible: status=%d\n",status);
		fprintf(ferr,"ERROR...  HeuSwaps solution not feasible: status=%d\n",status);
		//exit(status);
	}
	else if (((*zheu)>zheuf+LPrec)||((*zheu)<zheuf-LPrec))
	{
		printf("ERROR...  HeuSwaps value: zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		fprintf(ferr,"ERROR...  HeuSwaps value: zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		//exit(status);
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Compute an feasible solution using a Lagrangian Heuristic
//-----------------------------------------------------------------------------
//
int AgileOpt::LagrHeu(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheu, double *yheu)
{
	long i,j,k,h,j1,kk,k1;
	long kount;
	double ptot;
	double objtot;
	double objval;
	double zheut;
	double zheuf;
	double minur;
	double mina;
	int status,bolf;

	for (j=0; j<n; j++)
	{
		Flag[j] = 0;
	}

	// Setup the initial heuristic solution value  
	(*zheu) = 0.0;
	objval = 0.0;

	for (i=0; i<m; i++)
	{
		// LP local into the loop body
		CplexObj LP;

		// Set problem size
		LP.ncols = n+nYMap;
		LP.nrows = 1 + 2*nYMap;  // OR/AND constraints are added later
		LP.nz = (n)*(2+n) + (nYMap)*(1+n);  // We overestimate |Yj|=n
		if (cfg.Sentinel)
		{
			LP.ncols += 1;
			LP.nrows += 1;
			LP.nz += n+1;
		}
		LP.MallocCols();

		// Compute min urpen, in order to avoid negative price
		minur=0.0;
		mina=0.0;
		for (j=0; j<n; j++)
		{
			kk = i*n+j;
			if (minur > urpen[kk])
				minur = urpen[kk];
			if (mina > apen[kk])
				mina = apen[kk];
		}

		// Load Columns/Variables xj in the Cplex data structure
		k = 0;
		kount = 0;
		for (j=0; j<n; j++)
		{
			kk = i*n+j;
			//LP.obj[k] = (m-i)*urpen[k];
			LP.obj[k] = urpen[kk]-minur+10.;

			LP.matbeg[k] = kount;

			// Knapsack Constraint
			LP.matind[kount] = 0;
			// LP.matval[kount] = p[j]*run[j]; // or pr[j]
			LP.matval[kount] = pr[j];
			kount++;

			// Linking Constraints xj and yj
			for (h=0; h<nY[j]; h++)
			{
				j1 = Y[j][h];
				if ((YMap[j1]>=0)&&(j1!=j))
				{
					LP.matind[kount] = 1+YMap[j1];
					LP.matval[kount] = -1.;
					kount++;
				}
			}

			if (YMap[j]>=0)
			{
				LP.matind[kount] = 1+nYMap+YMap[j];
				LP.matval[kount] = -(double)(nY[j]-1);
				kount++;
			}

			// Sentinel Constraint
			if (cfg.Sentinel)
			{
				LP.matind[kount] = 1 + 2*nYMap;
				LP.matval[kount] = 1.0;
				kount++;
			}

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0;  // prima 1
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'B';

			LP.lb[k] = 0.0;
			if (Flag[j]==0)
				LP.ub[k] = 1.0;
			else
				LP.ub[k] = 0.0;
			k++;
		}

		// Load Columns/Variables yij in the Cplex data structure
		kk = 0;
		for (j=0; j<n; j++)
		{
			if (YMap[kk]<0)
			{
				kk++;
				continue;
			}

			k1 = i*n+j;
			if (nY[j]<2)
				LP.obj[k] = 0.0;
			else
				LP.obj[k] = apen[k1]-mina+10.;
				//LP.obj[k] = (m-i)*u[j]*a[j]/(nY[j]-1);
				// LP.obj[k] = a[j]/(nY[j]-1);

			LP.matbeg[k] = kount;

			LP.matind[kount] = 1+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matind[kount] = 1+nYMap+YMap[kk];
			LP.matval[kount] = +1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 0; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			LP.xctype[k] = 'C';
			//LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		if (cfg.Sentinel)
		{
			// Load Columns/Variables related to Sentinel constraints in the Cplex data structure
			kk = 0;

			LP.obj[k] = 0.0;

			LP.matbeg[k] = kount;

			LP.matind[kount] = m + n + 2*nYMap;
			LP.matval[kount] = -1.0;
			kount++;

			LP.matcnt[k] = kount - LP.matbeg[k];

			LP.indices[k] = k;
			LP.priority[k] = 1000; // Prima 0
			LP.direction[k] = CPX_BRANCH_GLOBAL;
			//LP.direction[k] = CPX_BRANCH_DOWN;
			LP.xctype[k] = 'I';

			LP.lb[k] = 0.0;
			LP.ub[k] = (double)n;
			k++;
			kk++;
		}

		LP.matbeg[k]=kount;

		// Load Constraints in the Cplex data structure
		kount=0;
		LP.rhs[kount]=pmax[i];
		LP.sense[kount]='L';
		kount++;

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		for (k=0; k<nYMap; k++)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='L';
			kount++;
		}

		if (cfg.Sentinel)
		{
			LP.rhs[kount]=0.0;
			LP.sense[kount]='E';
			kount++;
		}

		// Maximization problem
		LP.minmax=-1;

		// Load Problem
		LP.CopyLP();

		// Set MIP
		LP.SetMIP(cfg.TimeLimit);

		// Setup Custom Cutting Plane
		LP.CuttingPlaneHeu(n,m,ur,pr,pmax,nY,Y,nUOR,UOR,nUAND,UAND,Flag,FDepP,cfg.Pred,cfg.Cover,cfg.Lifting,cfg.KnapSol,cfg.MaxCuts,&Cuts,0);

		// Optimize
		LP.SolveMIP(&zheut,&Gap,&Nodes,&Cuts);

		//(*zheu) += zheut;

		// Read the solution
		ptot=0.0;
		objtot=0.0;
		if (iprinx>0) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);

		// Variables xj
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (LP.xt[j]>0.001)
			{
				xheu[k] = 1.0;
				(*zheu) += (m-i)*ur[j]*LP.xt[j];
				Flag[j] = 1;
				ptot += p[j]*run[j];
				objtot += LP.obj[j];
				if (iprinx>0) printf("%d ",j);
			}
			else
				xheu[k] = 0.0;
		}

		// Variables yj
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			if (YMap[j]<0)
			{
				yheu[k] = 0.0;
				continue;
			}

			kk = n + YMap[j];
			if (LP.xt[kk]>0.001)
			{
				yheu[k] = LP.xt[kk];
				(*zheu) += (m-i)*u[j]*(a[j]/((double)(nY[j]-1)))*LP.xt[kk];
				objtot += (m-i)*u[j]*(a[j]/((double)(nY[j]-1)))*LP.xt[kk];
			}
			else
				yheu[k] = 0.0;
		}

		objval += objtot;
		if (iprinx>0) printf("\n  ptot=%lf    ovjtot=%lf\n",ptot,objtot);
	}

	// Check feasibility
	status = Feasibility(xheu,yheu,&zheuf,&bolf);
	if (status)
	{
		printf("ERROR...  LagrHeuristic solution (MIP) not feasible: status=%d\n",status);
		fprintf(ferr,"ERROR...  LagrHeuristic solution (MIP) not feasible: status=%d\n",status);
		*zheu=-Inf;
		//exit(status);
	}
	else if (((*zheu)>zheuf+LPrec)||((*zheu)<zheuf-LPrec))
	{
		printf("ERROR...  LagrHeuristic value (MIP): zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		fprintf(ferr,"ERROR...  LagrHeuristic value (MIP): zheu=%lf  zheut=%lf\n",*zheu,zheuf);
		*zheu=-Inf;
		//exit(status);
	}

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zheu));

	return 0;
}


//-----------------------------------------------------------------------------
// Compute an feasible solution using a Lagrangian Heuristic
//-----------------------------------------------------------------------------
//
int AgileOpt::LagrHeuDP(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheub, double *yheub)
{
	long i,j,k,h,j1,jj,kk;
	long iter;
	long kount;
	long bol;
	double ptot;
	double objtot;
	double objval;
	double zheuf;
	double zheut;
	double zheub;
	double minur;
	double maxur;
	long jmax;
	long maxw;
	int status;
	int bolf;
	Knapsack KP;
	double *xheu;
	double *yheu;

	//(*zheu) = 0.0;
	(*zheu) = -Inf;
	xheu = new double [n*m];
	yheu = new double [n*m];

	maxw = 0;
	for (i=0; i<m; i++)
	{
		if (maxw<(int(pmax[i]*10+0.001)))
			maxw = int(pmax[i]*10+0.001);
	}

	// Allocate Knapsack data structures
	KP.Malloc(n,maxw);

	for (iter=0; iter<3; iter++)
	{
		for (j=0; j<n; j++)
		{
			Flag[j] = 0;
		}

		//maxw = 0;
		//for (i=0; i<m; i++)
		//{
		//	if (maxw<(int(pmax[i]*10+0.001)))
		//		maxw = int(pmax[i]*10+0.001);
		//}

		//// Allocate Knapsack data structures
		//KP.Malloc(n,maxw);

		// Setup the initial heuristic solution value  
		zheub = 0.0;
		objval = 0.0;

		for (i=0; i<m; i++)
		{
			// LP local into the loop body
			//CplexObj LP;

			// Compute min urpen, in order to avoid negative price
			minur=0.0;
			for (j=0; j<n; j++)
			{
				kk = i*n+j;
				if (minur > urpen[kk])
					minur = urpen[kk];
			}

retry:
			// Variables xij
			zheut = 0.0;
			KP.n = 0;
			KP.nf = 0;
			KP.zf = 0.0;
			KP.W = int(pmax[i]*10+0.001);
			for (j=0; j<n; j++) 
			{
				k = n*i+j;
				if (Flag[j]==0)
				{
					KP.w[KP.n] = int(pr[j]*10+0.001);
					KP.p[KP.n] = urpen[k]-minur+10.;
					KP.ind[KP.n] = j;
					KP.n++;
				}
				else if (Flag[j]==2)
				{
					KP.W -= int(pr[j]*10+0.001);
					KP.wf[KP.nf] = int(pr[j]*10+0.001);
					KP.zf += urpen[k]-minur+10.;
					KP.indf[KP.nf] = j;
					KP.nf++;
				}
			}

			if (KP.W < -Prec)
			{
				zheub = 0.0;
				break;
			}

			KP.SolveReMap();

			// Read the solution
			ptot=0.0;
			objtot=0.0;
			if (iprinx>1) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				xheu[k] = 0.0;
			}
			for (j1=0;j1<KP.nf;j1++)
			{
				j = KP.indf[j1];
				k = i*n+j;
				xheu[k] = 1.0;
				zheut += xheu[k] * (m-i)*u[j]*rcr[j];
				ptot += p[j]*run[j];
				objtot += xheu[k] * (m-i)*u[j]*rcr[j];
				if (iprinx>1) printf("%d ",j);
			}
			for (j1=0;j1<KP.n;j1++)
			{
				j = KP.ind[j1];
				k = i*n+j;
				xheu[k] = (double)KP.x[j1];
				if (xheu[k]>Prec)
				{
					zheut += xheu[k] * (m-i)*u[j]*rcr[j];
					ptot += p[j]*run[j];
					objtot += xheu[k] * (m-i)*u[j]*rcr[j];
					if (iprinx>1) printf("%d ",j);
				}
			}
			if (iprinx>1) printf("\n  ptot=%lf\n",ptot);

			// Feasibility check!
			// The relaxed precedence constraints must be checked
			for (j=0;j<n;j++)
			{
				k = i*n+j;
				if (xheu[k]<Prec)
					continue;

				goto ConAND;

ConOR:
				// OR Constraints
				if (nUOR[j]>0)
				{
					bol = 0;
					jmax = -1;
					maxur = -Inf;
					for (kk=0;kk<nUOR[j];kk++)
					{
						jj = UOR[j][kk];
						k = n*i+jj;
						if ((Flag[jj]==1)||(xheu[k]>0.999))
						{
							bol = 1;
							break;
						}
						if (maxur<(urpen[k]/pr[jj]))
						{
							jmax = jj;
							maxur = urpen[k]/pr[jj];
						}
					}
					if (bol==0)
					{
						if (iter==0)
							Flag[j] = -1;
						else 
							Flag[jmax] = +2;
						goto retry;
					}
				}

				goto Fine;

ConAND:
				// AND Constraints
				if (nUAND[j]>0)
				{
					jmax = -1;
					maxur = -Inf;
					bol = 0;
					for (kk=0;kk<nUAND[j];kk++)
					{
						jj = UAND[j][kk];
						k = n*i+jj;
						if ((Flag[jj]!=1)&&(xheu[k]<0.999))
						{
							if (iter==1)
								Flag[jj] = +2;
							else if (iter==2)
							{
								if (maxur<(urpen[k]/pr[jj]))
								{
									jmax = jj;
									maxur = urpen[k]/pr[jj];
								}
							}

							bol=1;
						}
						if (bol)
						{
							if (iter==0)
								Flag[j] = -1;
							else if (iter==2)
								Flag[jmax] = +2;
							goto retry;
						}
					}
				}
				goto ConOR;
Fine: ;
			}

			// Setta Flag
			for (j=0;j<n;j++)
			{
				if (Flag[j]<0)
					Flag[j] = 0;
				k = n*i+j;
				if (xheu[k]>0.999) 
					Flag[j] = 1;
			}

			ptot = 0.0;
			if (iprinx>0) printf("Sprint %d (pmax=%lf, zheut=%lf): ",i,pmax[i],zheut);
			zheut = 0.0;
			for (j=0;j<n;j++)
			{
				k = i*n+j;
				if (xheu[k]>Prec)
				{
					zheut += xheu[k] * (m-i)*u[j]*rcr[j];
					ptot += p[j]*run[j];
					if (iprinx>0) printf("%d ",j);
				}
			}

			// Variables yij
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				yheu[k] = 0.0;

				if ((nY[j]<2)||(xheu[k]<Prec))
					continue;

				//ptot = 0.0;
				for (j1=0; j1<nY[j]; j1++)
				{
					jj = Y[j][j1];
					if (jj==j)
						continue;
					kk = n*i+jj;
					yheu[k] += xheu[kk];
					//ptot += xheu[k];
				}

				zheut += yheu[k] * (m-i)*u[j]*a[j]/(nY[j]-1); 
			}

			zheub += zheut;
			if (iprinx>0) printf(" --> zheut0=%lf \n",zheut);

			// There are other items for the next sprints?
			//bol = 0;
			//for (j=0;j<n;j++)
			//{
			//	if (Flag[j]==0)
			//	{
			//		bol = 1;
			//		break;
			//	}
			//}
			//if (bol==0)
			//	break;
		}

		if ((*zheu)<zheub)
		{
			status = Feasibility(xheu,yheu,&zheuf,&bolf);
			if (status)
			{
				printf("ERROR...  Heuristic solution not feasible (SKIP): iter=%d  status=%d  zheu=%lf  \n",iter,status,(zheub));
				fprintf(ferr,"ERROR...  Heuristic solution not feasible (SKIP): iter=%d  status=%d  zheu=%lf  \n",iter,status,(zheub));
				//goto skip;
				//exit(status);
				continue;
			}
			else if ((zheub>zheuf+LPrec)||(zheub<zheuf-LPrec))
			{
				printf("ERROR...  Heuristic value (SKIP): iter=%d  zheu=%lf  zheut=%lf\n",iter,zheub,zheuf);
				fprintf(ferr,"ERROR...  Heuristic value (SKIP): iter=%d  zheu=%lf  zheut=%lf\n",iter,zheub,zheuf);
				//goto skip;
				//exit(status);
				continue;
			}

			(*zheu) = zheub; 
			for (j=0; j<n; j++)
				for (i=0; i<m; i++)
				{
					k = i*n+j;
					xheub[k] = xheu[k];
					yheub[k] = yheu[k];
				}
		}
	}

//skip:

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zheu));

	delete [] xheu;
	delete [] yheu;

	return 0;
}


//-----------------------------------------------------------------------------
// Compute an feasible solution using a Lagrangian Heuristic
//-----------------------------------------------------------------------------
//
int AgileOpt::LagrHeuDPBack(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheub, double *yheub)
{
	long i,j,k,h,j1,jj,kk;
	long iter,iter1;
	long kount;
	long bol;
	double ptot;
	double objtot;
	double objval;
	double zheuf;
	double zheut;
	double zheub;
	double minur;
	double maxur;
	long jmax;
	long maxw;
	int status;
	int bolf;
	Knapsack KP;
	double *xheu;
	double *yheu;
	double *corr;
	double *bcorr;
	double coef;

	//(*zheu) = 0.0;
	(*zheu) = -Inf;
	xheu = new double [n*m];
	yheu = new double [n*m];
	corr = new double [n];
	bcorr = new double [n];

	maxw = 0;
	for (i=0; i<m; i++)
	{
		if (maxw<(int(pmax[i]*10+0.001)))
			maxw = int(pmax[i]*10+0.001);
	}

	// Allocate Knapsack data structures
	KP.Malloc(n,maxw);

	for (iter=0; iter<5; iter++)
	{
		iter1 = 0;

		if (iter<3)
		{
			coef=1.0;
			for (j=0;j<n;j++)
				corr[j]=1.0;
		}
		else if (iter<5)
		{
			if (iter==3)
				coef=2.0;
			else
				coef=5.0;
			for (j=0;j<n;j++)
				corr[j]=1.0;
				//corr[j]=1.025;
		}
		else
		{
			coef=1.0;
			for (j=0;j<n;j++)
				corr[j]=1.025;
		}

retry4:
		iter1++;
		if (iter1>1000)
			break;
		bol=0;
		for (j=0;j<n;j++)
			if (corr[j]>+1.0e+20) 
			{
				bol=1;
				break;
			}
		if (bol)
			break;

		for (j=0; j<n; j++)
		{
			Flag[j] = 0;
		}

		// Setup the initial heuristic solution value  
		zheub = 0.0;
		objval = 0.0;

		for (i=0; i<m; i++)
		{
			// LP local into the loop body
			//CplexObj LP;

			// Compute min urpen, in order to avoid negative price
			minur=0.0;
			for (j=0; j<n; j++)
			{
				kk = i*n+j;
				if (minur > urpen[kk])
					minur = urpen[kk];
			}

retry:
			// Variables xij
			zheut = 0.0;
			KP.n = 0;
			KP.nf = 0;
			KP.zf = 0.0;
			KP.W = int(pmax[i]*10+0.001);
			for (j=0; j<n; j++) 
			{
				k = n*i+j;
				if (Flag[j]==0)
				{
					KP.w[KP.n] = int(pr[j]*10+0.001);
					KP.p[KP.n] = corr[j]*(urpen[k]-minur+10.);
					KP.ind[KP.n] = j;
					KP.n++;
				}
				else if (Flag[j]==2)
				{
					KP.W -= int(pr[j]*10+0.001);
					KP.wf[KP.nf] = int(pr[j]*10+0.001);
					KP.zf += urpen[k]-minur+10.;
					KP.indf[KP.nf] = j;
					KP.nf++;
				}
			}

			if (KP.W < -Prec)
			{
				zheub = 0.0;
				break;
			}

			KP.SolveReMap();

			// Read the solution
			ptot=0.0;
			objtot=0.0;
			if (iprinx>1) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				xheu[k] = 0.0;
			}
			for (j1=0;j1<KP.nf;j1++)
			{
				j = KP.indf[j1];
				k = i*n+j;
				xheu[k] = 1.0;
				zheut += xheu[k] * (m-i)*u[j]*rcr[j];
				ptot += p[j]*run[j];
				objtot += xheu[k] * (m-i)*u[j]*rcr[j];
				if (iprinx>1) printf("%d ",j);
			}
			for (j1=0;j1<KP.n;j1++)
			{
				j = KP.ind[j1];
				k = i*n+j;
				xheu[k] = (double)KP.x[j1];
				if (xheu[k]>Prec)
				{
					zheut += xheu[k] * (m-i)*u[j]*rcr[j];
					ptot += p[j]*run[j];
					objtot += xheu[k] * (m-i)*u[j]*rcr[j];
					if (iprinx>1) printf("%d ",j);
				}
			}
			if (iprinx>1) printf("\n  ptot=%lf\n",ptot);

			// Feasibility check!
			// The relaxed precedence constraints must be checked
			for (j=0;j<n;j++)
			{
				k = i*n+j;
				if (xheu[k]<Prec)
					continue;

				goto ConAND;

ConOR:
				// OR Constraints
				if (nUOR[j]>0)
				{
					bol = 0;
					jmax = -1;
					maxur = -Inf;
					for (kk=0;kk<nUOR[j];kk++)
					{
						jj = UOR[j][kk];
						k = n*i+jj;
						if ((Flag[jj]==1)||(xheu[k]>0.999))
						{
							bol = 1;
							break;
						}
						if (maxur<(urpen[k]/pr[jj]))
						{
							jmax = jj;
							maxur = urpen[k]/pr[jj];
						}
					}
					if (bol==0)
					{
						if (iter>=3)
						{
							for (kk=0;kk<nUOR[j];kk++)
							{
								jj = UOR[j][kk];
								k = n*i+jj;
								if ((Flag[jj]!=1)&&(xheu[k]<0.999))
								{
									if (iter<5)
										corr[jj] = coef*corr[jj];
									else
										corr[jj] = corr[jj]*corr[jj];
								}
							}
							goto retry4;
						}
						else if (iter==0)
							Flag[j] = -1;
						else 
							Flag[jmax] = +2;

						goto retry;
					}
				}

				goto Fine;

ConAND:
				// AND Constraints
				if (nUAND[j]>0)
				{
					jmax = -1;
					maxur = -Inf;
					bol = 0;
					for (kk=0;kk<nUAND[j];kk++)
					{
						jj = UAND[j][kk];
						k = n*i+jj;
						if ((Flag[jj]!=1)&&(xheu[k]<0.999))
						{
							if (iter>=3)
							{
								if (iter<5)
									corr[jj] = coef*corr[jj];
								else
									corr[jj] = corr[jj]*corr[jj];
							}
							else if (iter==1)
								Flag[jj] = +2;
							else if (iter==2)
							{
								if (maxur<(urpen[k]/pr[jj]))
								{
									jmax = jj;
									maxur = urpen[k]/pr[jj];
								}
							}

							bol=1;
						}
						if (bol)
						{
							if (iter>=3)
								goto retry4;
							else if (iter==0)
								Flag[j] = -1;
							else if (iter==2)
								Flag[jmax] = +2;

							goto retry;
						}
					}
				}
				goto ConOR;
Fine: ;
			}

			// Setta Flag
			for (j=0;j<n;j++)
			{
				if (Flag[j]<0)
					Flag[j] = 0;
				k = n*i+j;
				if (xheu[k]>0.999) 
					Flag[j] = 1;
			}

			ptot = 0.0;
			if (iprinx>0) printf("Sprint %d (pmax=%lf, zheut=%lf): ",i,pmax[i],zheut);
			zheut = 0.0;
			for (j=0;j<n;j++)
			{
				k = i*n+j;
				if (xheu[k]>Prec)
				{
					zheut += xheu[k] * (m-i)*u[j]*rcr[j];
					ptot += p[j]*run[j];
					if (iprinx>0) printf("%d ",j);
				}
			}

			// Variables yij
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				yheu[k] = 0.0;

				if ((nY[j]<2)||(xheu[k]<Prec))
					continue;

				//ptot = 0.0;
				for (j1=0; j1<nY[j]; j1++)
				{
					jj = Y[j][j1];
					if (jj==j)
						continue;
					kk = n*i+jj;
					yheu[k] += xheu[kk];
					//ptot += xheu[k];
				}

				zheut += yheu[k] * (m-i)*u[j]*a[j]/(nY[j]-1); 
			}

			zheub += zheut;
			if (iprinx>0) printf(" --> zheut0=%lf \n",zheut);

			// There are other items for the next sprints?
			//bol = 0;
			//for (j=0;j<n;j++)
			//{
			//	if (Flag[j]==0)
			//	{
			//		bol = 1;
			//		break;
			//	}
			//}
			//if (bol==0)
			//	break;
		}

		// Try swap!!!
		status = Feasibility(xheu,yheu,&zheuf,&bolf);
		if (status)
			continue;
		HeuSwaps(&zheub,xheu,yheu);

		if ((*zheu)<zheub)
		{
			status = Feasibility(xheu,yheu,&zheuf,&bolf);
			if (status)
			{
				printf("ERROR...  Heuristic solution not feasible (SKIP): iter=%d  status=%d  zheu=%lf  \n",iter,status,(zheub));
				fprintf(ferr,"ERROR...  Heuristic solution not feasible (SKIP): iter=%d  status=%d  zheu=%lf  \n",iter,status,(zheub));
				//goto skip;
				//exit(status);
				continue;
			}
			else if ((zheub>zheuf+LPrec)||(zheub<zheuf-LPrec))
			{
				printf("ERROR...  Heuristic value (SKIP): iter=%d  zheu=%lf  zheut=%lf\n",iter,zheub,zheuf);
				fprintf(ferr,"ERROR...  Heuristic value (SKIP): iter=%d  zheu=%lf  zheut=%lf\n",iter,zheub,zheuf);
				//goto skip;
				//exit(status);
				continue;
			}

			(*zheu) = zheub; 
			for (j=0; j<n; j++)
			{
				bcorr[j]=corr[j];
				for (i=0; i<m; i++)
				{
					k = i*n+j;
					xheub[k] = xheu[k];
					yheub[k] = yheu[k];
				}
			}
		}
	}

//skip:

	//status = Feasibility(xheub,yheub,&zheub,&bolf);
	//if (status)
	//	status=status;

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zheu));

	//for (i=0; i<m; i++)
	//	for (j=0; j<n; j++)
	//	{
	//		k=n*i+j;
	//		urpen[k]=bcorr[j]*urpen[k];
	//		if (urpen[k]<-1.0e+10)
	//			urpen[k]=-1.0e+10;
	//		if (urpen[k]>+1.0e+10) 
	//			urpen[k]=+1.0e+10;
	//	}

	delete [] xheu;
	delete [] yheu;
	delete [] corr;
	delete [] bcorr;

	return 0;
}


//-----------------------------------------------------------------------------
// Compute an feasible solution using a Lagrangian Heuristic
//-----------------------------------------------------------------------------
//
int AgileOpt::LagrHeuDP_Old(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheu)
{
	long i,j,k,h,j1,jj,kk;
	long iter;
	long kount;
	long bol;
	double ptot;
	double zheut;
	double zheub;
	double minur;
	long maxw;
	long attempt;
	long maxattempt;
	double *uraux;
	double aux;
	Knapsack KP;

	uraux = new double [n*m];

	(*zheu) = 0.0;

	for (iter=0; iter<2; iter++)
	{
		if (iter==0)
			maxattempt = 0;
		else
			maxattempt = 5;

		for (j=0; j<n; j++)
		{
			Flag[j] = 0;
		}

		maxw = 0;
		for (i=0; i<m; i++)
		{
			if (maxw<(int(pmax[i]*10+0.001)))
				maxw = int(pmax[i]*10+0.001);
		}

		// Allocate Knapsack data structures
		KP.Malloc(n,maxw);

		// Setup the initial heuristic solution value  
		zheub = 0.0;

		for (i=0; i<m; i++)
		{
			// LP local into the loop body
			//CplexObj LP;

			// Compute min urpen, in order to avoid negative price
			minur=0.0;
			for (j=0; j<n; j++)
			{
				kk = i*n+j;
				if (minur > urpen[kk])
					minur = urpen[kk];
			}
			for (j=0; j<n; j++)
			{
				kk = i*n+j;
				uraux[kk] = urpen[kk]-minur+1.;
			}

			attempt = 0;
retry:
			attempt++;

			// Variables xij
			zheut = 0.0;
			KP.n = 0;
			KP.W = int(pmax[i]*10+0.001);
			for (j=0; j<n; j++) 
			{
				k = n*i+j;
				if (Flag[j]==0)
				{
					KP.w[KP.n] = int(pr[j]*10+0.001);
					//KP.p[KP.n] = urpen[k]-minur+10.;
					KP.p[KP.n] = uraux[k];
					KP.ind[KP.n] = j;
					KP.n++;
				}
			}

			KP.SolveReMap();

			// Read the solution
			ptot=0.0;
			if (iprinx>0) printf("Sprint %d (pmax=%lf): ",i,pmax[i]);
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				xheu[k] = 0.0;
			}
			for (j1=0;j1<KP.n;j1++)
			{
				j = KP.ind[j1];
				k = i*n+j;
				xheu[k] = (double)KP.x[j1];
				if (xheu[k]>Prec)
				{
					zheut += xheu[k] * (m-i)*u[j]*rcr[j];
					ptot += p[j]*run[j];
					if (iprinx>0) printf("%d ",j);
				}
			}
			if (iprinx>0) printf("\n  ptot=%lf\n",ptot);

			// Feasibility check!
			// The relaxed precedence constraints must be checked
			for (j1=0;j1<KP.n;j1++)
			{
				j = KP.ind[j1];
				if (KP.x[j1]>0) 
				{
					// OR Constraints
					if (nUOR[j]>0)
					{
						bol = 0;
						for (kk=0;kk<nUOR[j];kk++)
						{
							jj = UOR[j][kk];
							k = n*i+jj;
							if ((Flag[jj]==1)||(xheu[k]>0.999))
							{
								//urpen[k] = urpen[n*i+j];
								if (iter==1)
								{
									aux = (uraux[n*i+j]/pr[j])*pr[jj]+1.;
									if (aux<uraux[k]+0.999)
										aux = uraux[k]*1.2;
									uraux[k] = aux;
								}
								else if (iter==2)
									uraux[k] = uraux[n*i+j];
								bol = 1;
								break;
							}
						}
						if (bol==0)
						{
							if ((iter==0)||(attempt>maxattempt))
								Flag[j] = -1;
							else if (attempt>0)
							{
								k = n*i+j;
								uraux[k] = uraux[k] / 1.1;
								//if (urpen[k]>Prec)
								//	urpen[k] = urpen[k] / 1.1;
								//else
								//	urpen[k] = urpen[k] * 1.1;
							}
							//if (minur > urpen[k])
							//	minur = urpen[k];
							goto retry;
						}
					}

					// AND Constraints
					if (nUAND[j]>0)
					{
						bol = 0;
						for (kk=0;kk<nUAND[j];kk++)
						{
							jj = UAND[j][kk];
							k = n*i+jj;
							if ((Flag[jj]!=1)&&(xheu[k]<0.999))
							{
								//urpen[k] = urpen[n*i+j];
								if (iter==1)
								{
									aux = (uraux[n*i+j]/pr[j])*pr[jj]+1.;
									if (aux<uraux[k]+0.999)
										aux = uraux[k]*1.2;
									uraux[k] = aux;
								}
								else if (iter==2)
									uraux[k] = uraux[n*i+j];
								bol=1;
							}
							if (bol)
							{
								if ((iter==0)||(attempt>maxattempt))
									Flag[j] = -1;
								else if (attempt>0)
								{
									k = n*i+j;
									uraux[k] = uraux[k] / 1.1;
									//if (urpen[k]>Prec)
									//	urpen[k] = urpen[k] / 1.1;
									//else
									//	urpen[k] = urpen[k] * 1.1;
								}
								//if (minur > urpen[k])
								//	minur = urpen[k];
								goto retry;
							}
						}
					}

				}
			}

			// Setta Flag
			for (j=0;j<n;j++)
			{
				if (Flag[j]<0)
					Flag[j] = 0;
				k = n*i+j;
				if (xheu[k]>0.999) 
					Flag[j] = 1;
			}

			// Variables yij
			for (j=0; j<n; j++)
			{
				if ((nY[j]<2)||(xheu[k]<Prec))
					continue;

				ptot = 0.0;
				for (j1=0; j1<nY[j]; j1++)
				{
					jj = Y[j][j1];
					if (jj==j)
						continue;
					k = n*i+jj;
					ptot += xheu[k];
				}

				zheut += ptot * (m-i)*u[j]*a[j]/(nY[j]-1); 
			}

			zheub += zheut;

		}

		if ((*zheu)<zheub)
			(*zheu) = zheub; 
	}

	if (iprinx>0) printf("\n  Zheu=%lf\n",(*zheu));

	delete [] uraux;

	return 0;
}


//-----------------------------------------------------------------------------
// Check if a solution is feasible and return its value
// (return 0 if feasible; return >0 if unfeasible)
//-----------------------------------------------------------------------------
//
int AgileOpt::Feasibility(double *xheu, double *yheu, double *zheu, int *bol)
{
	long i,j,k;
	long i1,j1,k1;
	long ii,jj,kk;
	double aux;

	//Setup initial values
	if (zheu!=NULL) (*zheu) = 0.0;
	(*bol) = 0;

	// Capacity Constraints
	for (i=0; i<m; i++)
	{
		aux = 0.0;
		for (j=0; j<n; j++)
		{
			k = i*n+j;
			aux += pr[j]*xheu[k];
		}
		if (aux>(pmax[i]+Prec))
		{
			*bol = 1;
			return (*bol);
		}
	}

	// Assignment Constraints
	for (j=0; j<n; j++)
	{
		aux = 0.0;
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			aux += xheu[k];
		}
		if ((aux<(1.-Prec))||(aux>(1.+Prec)))
		{
			*bol = 2;
			return (*bol);
		}
	}

	// OR Constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]<1)
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			if (xheu[k]<Prec)
				continue;

			aux = -xheu[k];
			for (jj=0; jj<nUOR[j]; jj++)
			{
				j1 = UOR[j][jj];
				for (i1=0; i1<=i; i1++)
				{
					k1 = i1*n+j1;
					aux += xheu[k1];
				}
			}
			if (aux<-Prec)
			{
				*bol = 3;
				return (*bol);
			}
		}
	}

	// AND Constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]<1)
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			if (xheu[k]<Prec)
				continue;

			aux = -((double)nUAND[j])*xheu[k];
			for (jj=0; jj<nUAND[j]; jj++)
			{
				j1 = UAND[j][jj];
				for (i1=0; i1<=i; i1++)
				{
					k1 = i1*n+j1;
					aux += xheu[k1];
				}
			}
			if (aux<-Prec)
			{
				*bol = 4;
				return (*bol);
			}
		}
	}

	// Y1 Constraints
	for (j=0; j<n; j++)
	{
		if (nY[j]<2)
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			aux = yheu[k];
			for (jj=0; jj<nY[j]; jj++)
			{
				j1 = Y[j][jj];
				if (j!=j1)
				{
					k1 = i*n+j1;
					aux -= xheu[k1];
				}
			}
			if (aux>Prec)
			{
				*bol = 5;
				return (*bol);
			}
		}
	}

	// Y2 Constraints
	for (j=0; j<n; j++)
	{
		if (nY[j]<2)
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			aux = yheu[k] - ((double)(nY[j]-1))*xheu[k];

			if (aux>Prec)
			{
				*bol = 6;
				return (*bol);
			}
		}
	}

	// Evaluate objective function
	if (zheu!=NULL)
	{
		(*zheu) = 0.0;
		for (i=0; i<m; i++)
			for (j=0; j<n; j++)
			{
				k = i*n+j;
				(*zheu) += (double)(m-i)*u[j] * rcr[j] * xheu[k];
				if (nY[j]>1)
					(*zheu) += (double)(m-i)*u[j] * (a[j]/((double)(nY[j]-1)))*yheu[k];
			}
	}

	return (*bol);

}

//-----------------------------------------------------------------------------
// Update penalties
//-----------------------------------------------------------------------------
//
int AgileOpt::UpdatePen(double alpha, double zheu, double zlr, double *xlr, double *Lambda, double *LambdaOR, double *LambdaAND, 
	                    double *subgr, double *subgrOR, double *subgrAND)
{
	long i,j,k;
	long i1,j1,k1;
	long jj;
	double step,den;
	
	den = 0.0;

	// Assignment constraints
	for (j=0; j<n; j++)
	{
		subgr[j] = -1.0;
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgr[j] += xlr[k];
		}
		den += subgr[j];
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrOR[k] = -xlr[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUOR[j]; jj++)
				{
					j1 = UOR[j][jj];
					k1 = i1*n+j1;
					subgrOR[k] += xlr[k1];
				}
			}
			den += subgrOR[k];
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrAND[k] = -(double)nUAND[j]*xlr[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUAND[j]; jj++)
				{
					j1 = UAND[j][jj];
					k1 = i1*n+j1;
					subgrAND[k] += xlr[k1];
				}
			}
			den += subgrAND[k];
		}
	}

	// Compute step
	step = alpha*(zlr*.99)/(den*den);
	//step = alpha*(zlr-zheu)/(den*den);

	// Assignment constraints
	for (j=0; j<n; j++)
	{
		Lambda[j] += step*subgr[j];
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaOR[k] += step*subgrOR[k];
			if (LambdaOR[k]>0.0)
				LambdaOR[k] = 0.0;
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaAND[k] += step*subgrAND[k];
			if (LambdaAND[k]>0.0)
				LambdaAND[k] = 0.0;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Update penalties
//-----------------------------------------------------------------------------
//
int AgileOpt::UpdatePenDP(double alpha, double zheu, double zlr, double *xlr, double *ylr, double *Lambda, double *LambdaOR, double *LambdaAND, 
	                      double *LambdaY1, double *LambdaY2, double *subgr, double *subgrOR, double *subgrAND, double *subgrY1, double *subgrY2)
{
	long i,j,k;
	long i1,j1,k1;
	long jj;
	double step,den;
	
	den = 0.0;

	// Assignment constraints
	for (j=0; j<n; j++)
	{
		subgr[j] = -1.0;
		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgr[j] += xlr[k];
		}
		den += (subgr[j]*subgr[j]);
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrOR[k] = -xlr[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUOR[j]; jj++)
				{
					j1 = UOR[j][jj];
					k1 = i1*n+j1;
					subgrOR[k] += xlr[k1];
				}
			}
			den += (subgrOR[k]*subgrOR[k]);
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrAND[k] = -(double)nUAND[j]*xlr[k];
			for (i1=0; i1<=i; i1++)
			{
				for (jj=0; jj<nUAND[j]; jj++)
				{
					j1 = UAND[j][jj];
					k1 = i1*n+j1;
					subgrAND[k] += xlr[k1];
				}
			}
			den += (subgrAND[k]*subgrAND[k]);
		}
	}

	// Constraints Y1
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrY1[k] = -ylr[k];
			for (jj=0; jj<nY[j]; jj++)
			{
				j1 = Y[j][jj];
				if (j1==j)
					continue;

				k1 = i*n+j1;
				subgrY1[k] += xlr[k1];
			}
			den += (subgrY1[k]*subgrY1[k]);
		}
	}

	// Constraints Y2
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			subgrY2[k] = -ylr[k];
			subgrY2[k] += (double)(nY[j]-1)*xlr[k];
			den += (subgrY1[k]*subgrY1[k]);
		}
	}

	// Compute step
//2012	step = alpha*(zlr*.99)/den;
	step = alpha*(0.01*zlr)/den;
	//step = alpha*(zlr-zheu)/den;

	// Assignment constraints
	for (j=0; j<n; j++)
	{
		Lambda[j] += step*subgr[j];
	}

	// OR Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUOR[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaOR[k] += step*subgrOR[k];
			if (LambdaOR[k]>0.0)
				LambdaOR[k] = 0.0;
		}
	}

	// AND Precedence constraints
	for (j=0; j<n; j++)
	{
		if (nUAND[j]==0) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaAND[k] += step*subgrAND[k];
			if (LambdaAND[k]>0.0)
				LambdaAND[k] = 0.0;
		}
	}

	// Constraints Y1
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaY1[k] += step*subgrY1[k];
			if (LambdaY1[k]>0.0)
				LambdaY1[k] = 0.0;
		}
	}

	// Constraints Y2
	for (j=0; j<n; j++)
	{
		if (nY[j]<2) 
			continue;

		for (i=0; i<m; i++)
		{
			k = i*n+j;
			LambdaY2[k] += step*subgrY2[k];
			if (LambdaY2[k]>0.0)
				LambdaY2[k] = 0.0;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
//  Reduce the instance size                                    
//-----------------------------------------------------------------------------
//
int AgileOpt::Reduction(void)
{
	Knapsack KP;
	long i,j,j1;
	long maxw;
	double dw,dwmin;

	maxw = 0;
	for (i=0; i<m; i++)
	{
		if (maxw<int(pmax[i]*10+0.001))
			maxw = int(pmax[i]*10+0.001);
	}

	// Setup Knapsack data structure
	KP.Malloc(n,maxw);

	// Try to reduce the Sprint capacity 
	for (i=0; i<m; i++)
	{
		KP.n = 0;
		KP.W = int(pmax[i]*10+0.001);

		for (j=0; j<n; j++) 
		{
			KP.w[KP.n] = int(pr[j]*10+0.001);
			KP.p[KP.n] = KP.w[KP.n];
			KP.ind[KP.n] = j;
			KP.n++;
		}

#ifdef ReMap
		KP.SolveReMap();
#else
		KP.Solve();
#endif

		pmax[i] = pmax[i] - (double)(KP.W-KP.zopt)/10.;
	}

	// Try to increase the story weight 
	for (j1=0; j1<n; j1++) 
	{
		dwmin = Inf;
		for (i=0; i<m; i++)
		{
			KP.n = 0;
			KP.W = int(pmax[i]*10+0.001) - int(pr[j1]*10+0.001);

			for (j=0; j<n; j++) 
			{
				if (j==j1) continue;

				KP.w[KP.n] = int(pr[j]*10+0.001);
				KP.p[KP.n] = KP.w[KP.n];
				KP.ind[KP.n] = j;
				KP.n++;
			}

#ifdef ReMap
			KP.SolveReMap();
#else
			KP.Solve();
#endif

			dw = (double)(KP.W-KP.zopt);
			if (dwmin>dw)
				dwmin = dw;
		}
		pr[j1] = pr[j1] + dwmin/10.;
	}

	return 0;
}