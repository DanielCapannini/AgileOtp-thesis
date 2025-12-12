//-----------------------------------------------------------------------------
//  File: CplexObj.cpp                                                       
//                                                                           
//  Project: Cplex solver interface                   
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 24.08.2011                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CplexObj.h"
//#include "Knapsack.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
CplexObj::CplexObj(void)
{
	int param;
	int ivalue;
	double rvalue;

	env = CPXopenCPLEX(&status);
	if (status)
	{
		printf("\nError... CPXopenCPLEX...\n");
		env = NULL;
		return;
	}

	lp = CPXcreateprob(env, &status, "lp");
	if (status)
	{
		printf("\nError... CPXcreateprob...\n");
		lp = NULL;
		return;
	}

	param = CPX_PARAM_DATACHECK;
	ivalue = 1;
	//ivalue=0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return;

	param = CPX_PARAM_SCRIND;
	ivalue = 1;
	// ivalue=0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return;

	param = CPX_PARAM_MIPDISPLAY;
	ivalue = 2;
	// ivalue=0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return;

	param = CPX_PARAM_SIMDISPLAY;
	//ivalue=2;
	ivalue = 0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return;

	minmax = 1;     // Default: minimize
	ncols = 0;      // Number of columns/variables
	nrows = 0;      // Number of rows/constraints
	nz = 0;         // Number of non-zero coefficients into the constraint matrix
	objsen = +1;    // Minimization/Maximization Flag (+1=Min; -1=Max)
	obj = NULL;     // Objective Function coefficients
	rhs = NULL;     // Right Hand Side
	sense = NULL;   // Constraint Senses
	matbeg = NULL;  // Pointers to the beginning of each columns of the constraint matrix
	matcnt = NULL;  // Cardinality of each columns of the constraint matrix
	matind = NULL;  // Row index associated to each coefficient of the constraint matrix
	matval = NULL;  // Value of each coefficient of the constraint matrix
	lb = NULL;      // Lower Bound value of each variable
	ub = NULL;      // Upper Bound value of each variable
	indices = NULL;
	priority = NULL;
	direction = NULL;
	xctype = NULL;
	rngval = NULL;

	nsos = 0;        // Number of SOS (Special Ordered Set) 
	nsosnz = 0;      // Total number of membres in all SOS added
	typesos = NULL;  // Array containing SOS type (CPX_TYPE_SOS1='1', CPX_TYPE_SOS2='2')
	sosbeg = NULL;   // Pointers to the beginning of each SOS
	sosind = NULL;   // Index associated to the SOS member
	soswt = NULL;    // Weight associated to the SOS member

	// Output
	xt = NULL;      // Variable Values

	// Data structure used to add cutting planes
	cutinfo.FlagLP = 0;
	cutinfo.env = env;
	cutinfo.lp = lp;
	cutinfo.lpkp = NULL;
	cutinfo.x = NULL;
	cutinfo.beg = NULL;
	cutinfo.ind = NULL;
	cutinfo.val = NULL;
	cutinfo.rhs = NULL;
	cutinfo.Flag2I = NULL;
	cutinfo.Flag3I = NULL;
	cutinfo.Tupla2I = NULL;
	cutinfo.Tupla3I = NULL;
	cutinfo.Flag2IRev = NULL;
	cutinfo.Flag3IRev = NULL;
	cutinfo.Tupla2IRev = NULL;
	cutinfo.Tupla3IRev = NULL;

}

//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
CplexObj::~CplexObj(void)
{
	long i, j;

	//status = CPXwriteprob(env,lp,"myprob.lp","LP");

	// Close the Cplex Enviroment
	if (cutinfo.lpkp != NULL)
	{
		status = CPXfreeprob(env, &(cutinfo.lpkp));
		if (status)
		{
			printf("\nError... CPXfreeprob...\n");
			return;
		}
	}

	status = CPXfreeprob(env, &lp);
	if (status)
	{
		printf("\nError... CPXfreeprob...\n");
		return;
	}

	status = CPXcloseCPLEX(&env);
	if (status)
	{
		printf("\nError... CPXcloseCPLEX...\n");
		return;
	}

	free(obj);
	free(rhs);
	free(sense);
	free(matbeg);
	free(matcnt);
	free(matind);
	free(matval);
	free(lb);
	free(ub);
	free(xt);
	free(indices);
	free(priority);
	free(direction);
	free(xctype);

	free(typesos);
	free(sosbeg);
	free(sosind);
	free(soswt);

	free(cutinfo.x);
	free(cutinfo.beg);
	free(cutinfo.ind);
	free(cutinfo.val);
	free(cutinfo.rhs);

	if (cutinfo.Tupla2I != NULL)
	{
		for (i = 0; i < cutinfo.n; i++)
			if (cutinfo.Tupla2I[i] != NULL)
				free(cutinfo.Tupla2I[i]);
		free(cutinfo.Tupla2I);
	}
	if (cutinfo.Tupla3I != NULL)
	{
		for (i = 0; i < cutinfo.n; i++)
			if (cutinfo.Tupla3I[i] != NULL)
				free(cutinfo.Tupla3I[i]);
		free(cutinfo.Tupla3I);
	}
	if (cutinfo.Tupla2IRev != NULL)
	{
		for (i = 0; i < cutinfo.n; i++)
			if (cutinfo.Tupla2IRev[i] != NULL)
				free(cutinfo.Tupla2IRev[i]);
		free(cutinfo.Tupla2IRev);
	}
	if (cutinfo.Tupla3IRev != NULL)
	{
		for (i = 0; i < cutinfo.n; i++)
			if (cutinfo.Tupla3IRev[i] != NULL)
				free(cutinfo.Tupla3IRev[i]);
		free(cutinfo.Tupla3IRev);
	}
	if (cutinfo.Flag2I != NULL) free(cutinfo.Flag2I);
	if (cutinfo.Flag3I != NULL) free(cutinfo.Flag3I);
	if (cutinfo.Flag2IRev != NULL) free(cutinfo.Flag2IRev);
	if (cutinfo.Flag3IRev != NULL) free(cutinfo.Flag3IRev);

}


//-----------------------------------------------------------------------------
//  MallocCols: Allocate the "column" data structure
//-----------------------------------------------------------------------------
//
int CplexObj::MallocCols(void)
{
	obj = (double*)malloc(ncols * sizeof(double));
	rhs = (double*)malloc(nrows * sizeof(double));
	sense = (char*)malloc(nrows * sizeof(char));
	matbeg = (int*)malloc((ncols + 4) * sizeof(int));
	matcnt = (int*)malloc((ncols + 4) * sizeof(int));
	matind = (int*)malloc(nz * sizeof(int));
	matval = (double*)malloc(nz * sizeof(double));
	lb = (double*)malloc(ncols * sizeof(double));
	ub = (double*)malloc(ncols * sizeof(double));
	xt = (double*)malloc(ncols * sizeof(double));
	indices = (int*)malloc((ncols + 4) * sizeof(int));
	priority = (int*)malloc((ncols + 4) * sizeof(int));
	direction = (int*)malloc((ncols + 4) * sizeof(int));
	xctype = (char*)malloc((ncols + 4) * sizeof(char));

	typesos = (char*)malloc(nsos * sizeof(char));
	sosbeg = (int*)malloc((nsos + 4) * sizeof(int));
	sosind = (int*)malloc(ncols * sizeof(int));
	soswt = (double*)malloc(ncols * sizeof(double));

	return 0;
}


//-----------------------------------------------------------------------------
//  CopyLP
//-----------------------------------------------------------------------------
//
int CplexObj::CopyLP(void)
{
	// Load the LP instance
	status = CPXcopylp(env, lp, ncols, nrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, rngval);

	// Set problem: minimization or maximization
	CPXchgobjsen(env, lp, minmax);

	return status;
}


//-----------------------------------------------------------------------------
//  Set a Lower Bound for the MIP
//-----------------------------------------------------------------------------
//
int CplexObj::SetLB(double LB)
{
	int Param;
	double DValue;

	// Set Parameter
	Param = CPX_PARAM_CUTLO;
	DValue = LB - 1.;
	status = CPXsetdblparam(env, Param, DValue);

	return 0;
}


//-----------------------------------------------------------------------------
//  SetMIP
//-----------------------------------------------------------------------------
//
int CplexObj::SetMIP(double TLim)
{
	int Param;
	int Value;
	double DValue;

	// Change the problem from LP to MIP
	status = CPXchgprobtype(env, lp, CPXPROB_MILP);

	// Charge the variable type (the default is "continuous")
	status = CPXchgctype(env, lp, ncols, indices, xctype);

	// Set Parameter
	Param = CPX_PARAM_TILIM;
	DValue = TLim;
	status = CPXsetdblparam(env, Param, DValue);

	Param = CPX_PARAM_EPAGAP;
	DValue = 1.0e-15;
	status = CPXsetdblparam(env, Param, DValue);

	// Load SOS 
	//status = CPXaddsos(env,lp,nsos,nsosnz,typesos,sosbeg,sosind,soswt,NULL);

	// Load Priority Order 
	//status = CPXcopyorder(env,lp,ncols,indices,priority,direction);

	// Set Parameter
	Param = CPX_PARAM_COVERS;
	Value = 0;
	status = CPXsetintparam(env, Param, Value);

	// Cuts Factor
	//Param = CPX_PARAM_CUTSFACTOR;
	//DValue = 2.0;
	//status = CPXsetdblparam (env,Param,DValue);

	return status;
}


//-----------------------------------------------------------------------------
//  Sentinel Constraint
//-----------------------------------------------------------------------------
//
int CplexObj::Sentinel(long n, long m)
{
	int i, j, k, k1;
	int* beg = NULL;
	int* ind = NULL;
	char* sense = NULL;
	double* rhs = NULL;
	double* val = NULL;
	double rr;

	// Build the Integer Sentinel Costraint
	ind = (int*)malloc(ncols * sizeof(int));		if (ind == NULL) goto ESCI;
	val = (double*)malloc(ncols * sizeof(double));	if (val == NULL) goto ESCI;
	beg = (int*)malloc((m + 4) * sizeof(int));			if (beg == NULL) goto ESCI;
	sense = (char*)malloc(m * sizeof(char));		if (sense == NULL) goto ESCI;
	rhs = (double*)malloc(m * sizeof(double));		if (rhs == NULL) goto ESCI;

	//srand(0);

	k = 0;
	for (i = 0; i < m; i++)
	{
		rhs[i] = 0.0;
		sense[i] = 'E';
		beg[i] = k;
		for (j = 0; j < n; j++)
		{
			k1 = i * n + j;
			//rr = (double)rand()/RAND_MAX;
			//rr = 1;
			//if ((xctype[k1]=='B')&&(rr<0.5))
			if (xctype[k1] == 'B')
			{
				ind[k] = indices[k1];
				val[k] = 1.0;
				k++;
			}
		}
		ind[k] = ncols;
		val[k] = -1.0;
		k++;
		indices[ncols] = ncols;
		priority[ncols] = 1000;
		//direction[ncols]=CPX_BRANCH_GLOBAL;
		direction[ncols] = CPX_BRANCH_DOWN;
		xctype[ncols] = 'I';
		ncols++;
	}
	beg[m] = k;

	// Load the Integer Sentinel Constraints
	status = CPXaddrows(env, lp, m, m, k, rhs, sense, beg, ind, val, NULL, NULL);

	goto TERMINATE;

ESCI:
	fprintf(stderr, "No memory for solution values.\n");

TERMINATE:
	free(beg);
	free(ind);
	free(val);
	free(sense);
	free(rhs);

	return status;
}


//-----------------------------------------------------------------------------
//  SolveMIP
//-----------------------------------------------------------------------------
//
int CplexObj::SolveMIP(double* Zopt, double* Gap, long* Nodes, long* Cuts)
{
	int ncuts;
	double RGap;

	//status = CPXchgprobtype(env,lp,CPXPROB_MILP);
	//status = CPXchgctype(env,lp,ncols,indices,xctype);

	//status = CPXwriteprob(env,lp,"myprob.lp","LP");

	status = CPXmipopt(env, lp);

	status = CPXgetstat(env, lp);
	//if ((status!=CPXMIP_OPTIMAL)&&(status!=CPXMIP_OPTIMAL_TOL)&&(status!=CPXMIP_TIME_LIM_FEAS)) // Cplex 11.2
	if ((status != CPXMIP_OPTIMAL) && (status != CPXMIP_OPTIMAL_TOL) && (status != CPXMIP_FEASIBLE) && (status != CPXMIP_TIME_LIM_FEAS)) // Cplex 12.2
	{
		printf("Non Ammissibile...\n", status);
		objval = +0.0;
		RGap = 1.0;
		goto TERMINATE;
	}

	// Read the solution
	status = CPXgetobjval(env, lp, &objval);
	status = CPXgetx(env, lp, xt, 0, ncols - 1);
	status = CPXgetmiprelgap(env, lp, &RGap);

TERMINATE:

	// Ouput
	*Zopt = objval;
	*Gap = RGap * 100.;
	*Nodes = CPXgetnodecnt(env, lp);
	status = CPXgetnumcuts(env, lp, CPX_CUT_USER, &ncuts);
	//status = CPXgetnumcuts(env,lp,CPX_CUT_COVER,&ncuts);
	//*Cuts = ncuts;

	return status;
}


//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int CplexObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax,
	long* nY, long** Y, long* nUOR, long** UOR, long* nUAND, long** UAND, long* FDepA, long* FDepP,
	int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long* Cuts)
{
	long cur_numcols;

	// Set parameters

	// Assure linear mappings between the presolved and original models 
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if (status)  goto TERMINATE;

	// Turn on traditional search for use with control callbacks 
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if (status)  goto TERMINATE;

	// Let MIP callbacks work on the original model
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if (status)  goto TERMINATE;

	// Let MIP do not reduce LP
	status = CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if (status)  goto TERMINATE;

	// Let MIP do not presolve MIP
	//status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
	//if ( status )  goto TERMINATE;

	*Cuts = 0;
	cur_numcols = CPXgetnumcols(env, lp);
	cutinfo.FlagLP = 0;
	cutinfo.lp = lp;
	cutinfo.lpkp = NULL;
	cutinfo.numcols = cur_numcols;
	cutinfo.Cuts = Cuts;
	cutinfo.numtoadd = 0;
	cutinfo.maxcuts = MaxCuts;
	cutinfo.n = n;
	cutinfo.m = m;
	cutinfo.u = u;
	cutinfo.p = p;
	cutinfo.pmax = pmax;
	cutinfo.nY = nY;
	cutinfo.Y = Y;
	cutinfo.nUOR = nUOR;
	cutinfo.UOR = UOR;
	cutinfo.nUAND = nUAND;
	cutinfo.UAND = UAND;
	cutinfo.FDepA = FDepA;
	cutinfo.FDepP = FDepP;
	cutinfo.Pred = Pred;
	cutinfo.Cover = Cover;
	cutinfo.Lifting = Lifting;
	cutinfo.KnapSol = KnapSol;

	cutinfo.x = (double*)malloc(cur_numcols * sizeof(double));
	if (cutinfo.x == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.ind = (int*)malloc(cur_numcols * sizeof(int));
	if (cutinfo.ind == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.val = (double*)malloc(cur_numcols * sizeof(double));
	if (cutinfo.val == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.beg = (int*)malloc(2 * sizeof(int));
	if (cutinfo.beg == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.rhs = (double*)malloc(2 * sizeof(double));
	if (cutinfo.rhs == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.Flag2I = (int*)malloc(n * sizeof(int));
	if (cutinfo.Flag2I == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.Tupla2I = (struct tupla**)malloc(n * sizeof(struct tupla*));
	if (cutinfo.Tupla2I == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	CalcolaFlag2I(&cutinfo);
	cutinfo.Flag3I = (int*)malloc(n * sizeof(int));
	if (cutinfo.Flag3I == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.Tupla3I = (struct tupla**)malloc(n * sizeof(struct tupla*));
	if (cutinfo.Tupla3I == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	CalcolaFlag3I(&cutinfo);

	cutinfo.Flag2IRev = (int*)malloc(n * sizeof(int));
	if (cutinfo.Flag2IRev == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.Tupla2IRev = (struct tupla**)malloc(n * sizeof(struct tupla*));
	if (cutinfo.Tupla2IRev == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	CalcolaFlag2IRev(&cutinfo);
	cutinfo.Flag3IRev = (int*)malloc(n * sizeof(int));
	if (cutinfo.Flag3IRev == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.Tupla3IRev = (struct tupla**)malloc(n * sizeof(struct tupla*));
	if (cutinfo.Tupla3IRev == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	CalcolaFlag3IRev(&cutinfo);

	// Setup the data structure in Cplex for the Knapsack 
	if (KnapSol == 0)
		SetupKPDP(&cutinfo);
	else
		SetupKPMIP(&cutinfo);

	// Set up to use MIP callback 
	//status = CPXsetcutcallbackfunc (env, mycutcallback, &cutinfo);
	status = CPXsetusercutcallbackfunc(env, myNewcutcallback, &cutinfo);
	if (status)  goto TERMINATE;

TERMINATE:

	return 0;
}


//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int CplexObj::CuttingPlaneHeu(long n, long m, double* u, double* p, double* pmax,
	long* nY, long** Y, long* nUOR, long** UOR, long* nUAND, long** UAND, long* FDepA, long* FDepP,
	int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long* Cuts, long fact)
{
	long cur_numcols;
	int param;
	int ivalue;

	// Set parameters

	param = CPX_PARAM_SCRIND;
	//ivalue=1;
	ivalue = 0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return 1;

	param = CPX_PARAM_MIPDISPLAY;
	//ivalue=2;
	ivalue = 0;
	status = CPXsetintparam(env, param, ivalue);
	if (status > 0) return 1;

	// Check if cutting plane must be applied
	if (fact > 0)
		goto TERMINATE;

	// Assure linear mappings between the presolved and original models 
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if (status)  goto TERMINATE;

	// Turn on traditional search for use with control callbacks 
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if (status)  goto TERMINATE;

	// Let MIP callbacks work on the original model
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if (status)  goto TERMINATE;

	// Let MIP do not reduce LP
	status = CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if (status)  goto TERMINATE;

	// Let MIP do not presolve MIP
	//status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
	//if ( status )  goto TERMINATE;

	*Cuts = 0;
	cur_numcols = CPXgetnumcols(env, lp);
	cutinfo.FlagLP = 0;
	cutinfo.lp = lp;
	cutinfo.lpkp = NULL;
	cutinfo.numcols = cur_numcols;
	cutinfo.Cuts = Cuts;
	cutinfo.numtoadd = 0;
	cutinfo.maxcuts = MaxCuts;
	cutinfo.n = n;
	cutinfo.m = m;
	cutinfo.u = u;
	cutinfo.p = p;
	cutinfo.pmax = pmax;
	cutinfo.nY = nY;
	cutinfo.Y = Y;
	cutinfo.nUOR = nUOR;
	cutinfo.UOR = UOR;
	cutinfo.nUAND = nUAND;
	cutinfo.UAND = UAND;
	cutinfo.FDepA = FDepA;
	cutinfo.FDepP = FDepP;
	cutinfo.Pred = Pred;
	cutinfo.Cover = Cover;
	cutinfo.Lifting = Lifting;
	cutinfo.KnapSol = KnapSol;

	cutinfo.x = (double*)malloc(cur_numcols * sizeof(double));
	if (cutinfo.x == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.ind = (int*)malloc(cur_numcols * sizeof(int));
	if (cutinfo.ind == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.val = (double*)malloc(cur_numcols * sizeof(double));
	if (cutinfo.val == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.beg = (int*)malloc(2 * sizeof(int));
	if (cutinfo.beg == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.rhs = (double*)malloc(2 * sizeof(double));
	if (cutinfo.rhs == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}

	// Set up to use MIP callback 
	//status = CPXsetcutcallbackfunc (env, mycutcallback, &cutinfo);
	status = CPXsetusercutcallbackfunc(env, myHeucutcallback, &cutinfo);
	if (status)  goto TERMINATE;

TERMINATE:

	return 0;
}



int CplexObj::CuttingPlaneHeu_new(
	long n, long m,
	double* u, double* p, double* pmax,
	long* nY, long** Y,
	long* nUOR, long** UOR,
	long* nUAND, long** UAND,
	long* FDepA, long* FDepP,
	int Pred, int Cover, int Lifting, int KnapSol,
	int MaxCuts, long* Cuts, long fact)
{
	int status = 0;
	long cur_numcols = CPXgetnumcols(env, lp);

	*Cuts = 0;

	if (fact > 0) {
		return 0;
	}

	// Disabilita solo il presolve che modifica la dimensione del modello
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if (status) return 1;

	status = CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if (status) return 1;

	// Per sicurezza mantieni le dimensioni degli indici (serve per callback)
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if (status) return 1;

	cutinfo.lp = lp;
	cutinfo.numcols = cur_numcols;
	cutinfo.Cuts = Cuts;
	cutinfo.maxcuts = MaxCuts;

	cutinfo.FDepA = FDepA;
	cutinfo.FDepP = FDepP;

	cutinfo.u = u;
	cutinfo.p = p;
	cutinfo.pmax = pmax;
	cutinfo.nY = nY;
	cutinfo.Y = Y;

	cutinfo.nUOR = nUOR;
	cutinfo.UOR = UOR;
	cutinfo.nUAND = nUAND;
	cutinfo.UAND = UAND;

	cutinfo.Pred = Pred;
	cutinfo.Cover = Cover;
	cutinfo.Lifting = Lifting;
	cutinfo.KnapSol = KnapSol;

	cutinfo.x = (double*)malloc(cur_numcols * sizeof(double));
	cutinfo.ind = (int*)malloc(cur_numcols * sizeof(int));
	cutinfo.val = (double*)malloc(cur_numcols * sizeof(double));
	cutinfo.beg = (int*)malloc(sizeof(int));
	cutinfo.rhs = (double*)malloc(sizeof(double));
	cutinfo.last_cut_checksum = 0UL;
	cutinfo.last_activity = 0.0;
	cutinfo.last_rhs = 0.0;

	cutinfo.MaxCutsPerNode = 5;          // valore consigliato
	cutinfo.cuts_added_this_node = 0;

	if (!cutinfo.x || !cutinfo.ind || !cutinfo.val) {
		printf("Memory error in cutting plane data.\n");
		return 1;
	}

	status = CPXsetusercutcallbackfunc(env, myHeucutcallback_new, &cutinfo);
	if (status) return 1;

	return 0;
}

//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
static int CPXPUBLIC mycutcallback(CPXCENVptr env,
	void* cbdata,
	int        wherefrom,
	void* cbhandle,
	int* useraction_p)
{
	int status = 0;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;

	*useraction_p = CPX_CALLBACK_DEFAULT;

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols - 1);
	if (status)
	{
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}

	// OR precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUOR[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = i * n + j;
			if (x[k] < 0.001) continue;

			cutvio = -x[k];
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
				{
					k = i1 * n + cutinfo->UOR[j][j1];
					cutvio += x[k];
				}
			if (cutvio < -0.001)
			{
				k1 = 0;
				rhs[0] = 0.0;
				beg[0] = k1;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
					{
						k = i1 * n + cutinfo->UOR[j][j1];
						ind[k1] = k;
						val[k1] = 1.0;
						k1++;
					}
				k = i * n + j;
				ind[k1] = k;
				val[k1] = -1.0;
				k1++;
				beg[1] = k1;

				goto ADDCUT;
			}
		}
	}

	// AND precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUAND[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = i * n + j;
			if (x[k] < 0.001) continue;

			cutvio = -(double)cutinfo->nUAND[j] * x[k];
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
				{
					k = i1 * n + cutinfo->UAND[j][j1];
					cutvio += x[k];
				}
			if (cutvio < -0.001)
			{
				k1 = 0;
				rhs[0] = 0.0;
				beg[0] = k1;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
					{
						k = i1 * n + cutinfo->UAND[j][j1];
						ind[k1] = k;
						val[k1] = 1.0;
						k1++;
					}
				k = i * n + j;
				ind[k1] = k;
				val[k1] = -(double)cutinfo->nUAND[j];
				k1++;
				beg[1] = k1;

				goto ADDCUT;
			}
		}
	}

ADDCUT:
	// Use a cut violation tolerance of 0.01
	if (cutvio < -0.001)
	{
		// status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
		status = CPXcutcallbackadd(env, cbdata, wherefrom, k1, rhs[0], 'G', ind, val, 0);
		if (status)
		{
			fprintf(stderr, "Failed to add cut.\n");
			goto TERMINATE;
		}
	}

	// Tell CPLEX that cuts have been created
	*useraction_p = CPX_CALLBACK_SET;

TERMINATE:

	return (status);

}


//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
static int CPXPUBLIC myNewcutcallback(CPXCENVptr env,
	void* cbdata,
	int        wherefrom,
	void* cbhandle,
	int* useraction_p)
{
	long i;
	int status = 0;
	int BolAdd;
	char sense;
	int purgeable;
	int param;
	int ivalue;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Are we solving a subproblem?
	if (cutinfo->FlagLP)
		return (status);

	// Do we have reached the maximum number of cuts?  
	if (cutinfo->numtoadd > cutinfo->maxcuts)
		return (status);

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, cutinfo->x, 0, cutinfo->numcols - 1);
	if (status)
	{
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}

	//printf("\nSolution:\n");
	//for (i=0; i<cutinfo->numcols; i++)
	//	if (cutinfo->x[i]>0.00001)
	//		printf("x(%d)=%lf\n",i,cutinfo->x[i]);

	BolAdd = 0;

	BolAdd = OR_AND_Inequalities(cutinfo);
	if (BolAdd)
	{
		sense = 'G';
		purgeable = 0;
		goto AddCut;
	}

	if (cutinfo->Pred)
	{
		BolAdd = Precedence_Inequalities(cutinfo);
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
		BolAdd = Precedence_Inequalities_2Item(cutinfo);
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
		//BolAdd = Precedence_Inequalities_2ItemRev(cutinfo);
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
		BolAdd = Precedence_Inequalities_3Item(cutinfo);
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
		//BolAdd = Precedence_Inequalities_3ItemRev(cutinfo);
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
	}

	if (cutinfo->Cover)
	{
		if (cutinfo->KnapSol == 0)
			BolAdd = Cover_Inequalities(cutinfo);
		else
		{
			cutinfo->FlagLP = 1;

			param = CPX_PARAM_SCRIND;ivalue = 0;
			status = CPXsetintparam(cutinfo->env, param, ivalue);
			if (status > 0) return 1;

			param = CPX_PARAM_MIPDISPLAY;ivalue = 0;
			status = CPXsetintparam(cutinfo->env, param, ivalue);
			if (status > 0) return 1;

			BolAdd = Cover_InequalitiesCplex(cutinfo);

			param = CPX_PARAM_SCRIND;ivalue = 1;
			status = CPXsetintparam(cutinfo->env, param, ivalue);
			if (status > 0) return 1;

			param = CPX_PARAM_MIPDISPLAY;ivalue = 2;
			status = CPXsetintparam(cutinfo->env, param, ivalue);
			if (status > 0) return 1;

			cutinfo->FlagLP = 0;
		}
		if (BolAdd)
		{
			sense = 'L';
			purgeable = 1;
			goto AddCut;
		}
	}
AddCut:

	// Use a cut violation tolerance of 0.01
	if (BolAdd == 1)
	{
		cutinfo->numtoadd++;
		cutinfo->Cuts[0]++;
		if ((cutinfo->numtoadd % 1000) == 0) printf("\n  User Cuts added: %d\n\n", cutinfo->numtoadd);

		//if ((cutinfo->numtoadd%5000)==0) 
		//	printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);

		// status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
		status = CPXcutcallbackadd(env, cbdata, wherefrom, cutinfo->nz, cutinfo->rhs[0], sense, cutinfo->ind, cutinfo->val, purgeable);
		if (status)
		{
			fprintf(stderr, "Failed to add cut.\n");
			goto TERMINATE;
		}

		// Tell CPLEX that cuts have been created
		*useraction_p = CPX_CALLBACK_SET;

	}

TERMINATE:

	return (status);

}


//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
static int CPXPUBLIC myHeucutcallback(CPXCENVptr env,
	void* cbdata,
	int        wherefrom,
	void* cbhandle,
	int* useraction_p)
{
	long i;
	int status = 0;
	int BolAdd;
	char sense;
	int purgeable;
	int param;
	int ivalue;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	*useraction_p = CPX_CALLBACK_DEFAULT;

	// We are solving a subproblem?
	if (cutinfo->FlagLP)
		return (status);

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, cutinfo->x, 0, cutinfo->numcols - 1);
	if (status)
	{
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}

	//printf("\nSolution:\n");
	//for (i=0; i<cutinfo->numcols; i++)
	//	if (cutinfo->x[i]>0.00001)
	//		printf("x(%d)=%lf\n",i,cutinfo->x[i]);

	BolAdd = 0;

	BolAdd = OR_AND_InequalitiesHeu(cutinfo);
	if (BolAdd)
	{
		sense = 'G';
		purgeable = 0;
		goto AddCut;
	}

AddCut:

	// Use a cut violation tolerance of 0.01
	if (BolAdd == 1)
	{
		cutinfo->numtoadd++;
		cutinfo->Cuts[0]++;
		if ((cutinfo->numtoadd % 1000) == 0) printf("\n  User Cuts added: %d\n\n", cutinfo->numtoadd);

		// status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
		status = CPXcutcallbackadd(env, cbdata, wherefrom, cutinfo->nz, cutinfo->rhs[0], sense, cutinfo->ind, cutinfo->val, purgeable);
		if (status)
		{
			fprintf(stderr, "Failed to add cut.\n");
			goto TERMINATE;
		}

		// Tell CPLEX that cuts have been created
		*useraction_p = CPX_CALLBACK_SET;

	}

TERMINATE:

	return (status);

}

static int CPXPUBLIC myHeucutcallback_new(
	CPXCENVptr env,
	void* cbdata,
	int wherefrom,
	void* cbhandle,
	int* useraction_p)
{
	int status = 0;
	*useraction_p = CPX_CALLBACK_DEFAULT;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	// Non aggiungere tagli durante LP puro
	if (cutinfo->FlagLP)
		return 0;

	// Limita ai primi livelli dell'albero: tagli più efficaci
	int depth = 0;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
	if (depth > 3)
		return 0;

	// Recupera solo le colonne effettivamente usate nei tagli
	status = CPXgetcallbacknodex(env, cbdata, wherefrom, cutinfo->x, 0, cutinfo->numcols - 1);
	if (status) return status;

	// Trova più tagli e li aggiunge tutti
	int numCutsAdded = 0;
	printf("\n mi sono fermato qui 1");

	while (1) {
		int violated = OR_AND_InequalitiesHeu(cutinfo);
		printf("\n mi sono fermato qui 2");

		if (!violated)
			break;

		status = CPXcutcallbackadd(env, cbdata, wherefrom,
			cutinfo->nz,
			cutinfo->rhs[0],
			'G',
			cutinfo->ind,
			cutinfo->val,
			CPX_USECUT_FILTER); // permette a CPLEX di eliminarli

		if (status) return status;
		printf("\n mi sono fermato qui 3");
		numCutsAdded++;
	}
	printf("\n mi sono fermato qui 4");
	if (numCutsAdded > 0)
		*useraction_p = CPX_CALLBACK_SET;

	return 0;
}


int OR_AND_Inequalities_new(CUTINFOptr cutinfo)
{
	const double EPS_VIOL = 1e-6;   // tolleranza per considerare una violazione reale
	int status = 0;

	int      numcols = cutinfo->numcols;
	double* x = cutinfo->x;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;

	long i, j, i1, j1;
	int k;
	int k1;

	// piccoli helper
	auto make_checksum = [&](int* ids, int nn)->unsigned long {
		unsigned long h = 1469598103934665603UL; // FNV-1a 64bit start
		for (int t = 0;t < nn;t++) {
			h ^= (unsigned long)ids[t];
			h *= 1099511628211UL;
		}
		return h;
		};

	// OR precedence
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUOR[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = (int)(i * n + j);
			if (x[k] < 1e-9) continue; // non attivo

			// calcola lhs = sum_{i1<=i, u in UOR[j]} x[i1,u]  e rhs = x[i,j]
			double lhs = 0.0;
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
				{
					int col = (int)(i1 * n + cutinfo->UOR[j][j1]);
					lhs += x[col];
				}
			double right = x[k];

			double viol = (lhs - right); // l'hai costruito come lhs - rhs >= 0 aspettato
			if (viol > EPS_VIOL) {
				// costruisci cut: sum_{i1<=i,u in UOR[j]} x[i1,u] - x[i,j] <= 0
				k1 = 0;
				rhs[0] = 0.0;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
					{
						int col = (int)(i1 * n + cutinfo->UOR[j][j1]);
						ind[k1] = col;
						val[k1] = 1.0;
						k1++;
					}
				ind[k1] = k;
				val[k1] = -1.0;
				k1++;

				// semplice controllo anti-duplicato rapido: checksum delle colonne
				unsigned long checksum = make_checksum(ind, k1);
				if (cutinfo->last_cut_checksum == checksum) {
					// già aggiunto negli ultimi cut processati: skip
					return 0;
				}
				cutinfo->last_cut_checksum = checksum;

				cutinfo->nz = k1;
				cutinfo->beg[0] = 0;
				cutinfo->beg[1] = k1;

				// opzionale: salva anche activity per controllo esterno/debug
				cutinfo->last_activity = lhs;
				cutinfo->last_rhs = right;

				return 1;
			}
		}
	}

	// AND precedence
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUAND[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = (int)(i * n + j);
			if (x[k] < 1e-9) continue;

			// lhs = sum_{i1<=i, u in UAND[j]} x[i1,u]
			double lhs = 0.0;
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
				{
					int col = (int)(i1 * n + cutinfo->UAND[j][j1]);
					lhs += x[col];
				}
			double right = (double)cutinfo->nUAND[j] * x[k];

			double viol = (lhs - right); // must be >= 0
			if (viol > EPS_VIOL) {
				k1 = 0;
				rhs[0] = 0.0;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
					{
						int col = (int)(i1 * n + cutinfo->UAND[j][j1]);
						ind[k1] = col;
						val[k1] = 1.0;
						k1++;
					}
				ind[k1] = k;
				val[k1] = -(double)cutinfo->nUAND[j];
				k1++;

				unsigned long checksum = make_checksum(ind, k1);
				if (cutinfo->last_cut_checksum == checksum) {
					return 0;
				}
				cutinfo->last_cut_checksum = checksum;

				cutinfo->nz = k1;
				cutinfo->beg[0] = 0;
				cutinfo->beg[1] = k1;

				cutinfo->last_activity = lhs;
				cutinfo->last_rhs = right;

				return 1;
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for OR-AND constraints
//-----------------------------------------------------------------------------
//
int OR_AND_Inequalities(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;

	// OR precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUOR[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = i * n + j;
			if (x[k] < 0.001) continue;

			cutvio = -x[k];
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
				{
					k = i1 * n + cutinfo->UOR[j][j1];
					cutvio += x[k];
				}
			if (cutvio < -0.001)
			{
				k1 = 0;
				rhs[0] = 0.0;
				beg[0] = k1;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
					{
						k = i1 * n + cutinfo->UOR[j][j1];
						ind[k1] = k;
						val[k1] = 1.0;
						k1++;
					}
				k = i * n + j;
				ind[k1] = k;
				val[k1] = -1.0;
				k1++;
				beg[1] = k1;

				cutinfo->nz = k1;

				return 1;
			}
		}
	}

	// AND precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUAND[j] == 0) continue;

		for (i = 0; i < m; i++)
		{
			k = i * n + j;
			if (x[k] < 0.001) continue;

			cutvio = -(double)cutinfo->nUAND[j] * x[k];
			for (i1 = 0; i1 <= i; i1++)
				for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
				{
					k = i1 * n + cutinfo->UAND[j][j1];
					cutvio += x[k];
				}
			if (cutvio < -0.001)
			{
				k1 = 0;
				rhs[0] = 0.0;
				beg[0] = k1;
				for (i1 = 0; i1 <= i; i1++)
					for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
					{
						k = i1 * n + cutinfo->UAND[j][j1];
						ind[k1] = k;
						val[k1] = 1.0;
						k1++;
					}
				k = i * n + j;
				ind[k1] = k;
				val[k1] = -(double)cutinfo->nUAND[j];
				k1++;
				beg[1] = k1;

				cutinfo->nz = k1;

				return 1;
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for OR-AND constraints
//-----------------------------------------------------------------------------
//
int OR_AND_InequalitiesHeu(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;

	// OR precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUOR[j] == 0) continue;

		if (x[j] < 0.001) continue;

		cutvio = -x[j];
		for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
		{
			k = cutinfo->UOR[j][j1];
			if (cutinfo->FDepA[k] > 0)
				cutvio += 1.0;
			else
				cutvio += x[k];
		}
		if (cutvio < -0.001)
		{
			k1 = 0;
			rhs[0] = 0.0;
			beg[0] = k1;
			for (j1 = 0; j1 < cutinfo->nUOR[j]; j1++)
			{
				k = cutinfo->UOR[j][j1];
				if (cutinfo->FDepA[k] == 0)
				{
					ind[k1] = k;
					val[k1] = 1.0;
					k1++;
				}
			}
			k = j;
			ind[k1] = k;
			val[k1] = -1.0;
			k1++;
			beg[1] = k1;

			cutinfo->nz = k1;

			return 1;
		}
	}

	// AND precedence
	cutvio = 0.0;
	for (j = 0; j < n; j++)
	{
		if (cutinfo->nUAND[j] == 0) continue;

		if (x[j] < 0.001) continue;

		cutvio = -(double)cutinfo->nUAND[j] * x[j];
		for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
		{
			k = cutinfo->UAND[j][j1];
			if (cutinfo->FDepA[k] > 0)
				cutvio += 1.0;
			else
				cutvio += x[k];
		}
		if (cutvio < -0.001)
		{
			k1 = 0;
			rhs[0] = 0.0;
			beg[0] = k1;
			for (j1 = 0; j1 < cutinfo->nUAND[j]; j1++)
			{
				k = cutinfo->UAND[j][j1];
				if (cutinfo->FDepA[k] == 0)
				{
					ind[k1] = k;
					val[k1] = 1.0;
					k1++;
				}
			}
			k = j;
			ind[k1] = k;
			val[k1] = -(double)cutinfo->nUAND[j];
			k1++;
			beg[1] = k1;

			cutinfo->nz = k1;

			return 1;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1, kk;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		if ((cutinfo->FDepA[j1] > 0) || (cutinfo->nY[j1] > 1))
			//if ((cutinfo->nY[j1]>1))
			continue;

		for (j = 0; j < n; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			if (fabs(cutinfo->p[j1] - cutinfo->p[j]) > 0.001)
				continue;

			if (cutinfo->u[j1] < cutinfo->u[j] - 0.001)
				continue;

			if ((fabs(cutinfo->u[j1] - cutinfo->u[j]) < 0.001) && (j1 > j))
				continue;

			for (i1 = 0; i1 < m; i1++)
			{
				k1 = i1 * n + j1;
				if (x[k1] < 0.001)
					continue;

				cutvio = 0.0;
				for (i = 0; i < i1 - 1; i++)
				{
					k = i * n + j;
					cutvio += x[k];
				}
				if ((cutvio + x[k1]) > 1.001)
				{
					//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
					//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
					//for (i1=0; i1<KP.n; i1++)
					//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

					kk = 0;
					rhs[0] = 1.0;
					beg[0] = kk;

					ind[kk] = k1;
					val[kk] = 1.0;
					kk++;

					for (i = 0; i < i1 - 1; i++)
					{
						k = i * n + j;
						ind[kk] = k;
						val[kk] = 1.0;
						kk++;
					}

					beg[1] = kk;
					cutinfo->nz = kk;

					return 1;
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_2Item(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, h, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag2I[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (h = 0; h < cutinfo->Flag2I[j1]; h++)
		{
			j = cutinfo->Tupla2I[j1][h].j1;
			j0 = cutinfo->Tupla2I[j1][h].j2;

			//if ((j==j1)||(cutinfo->FDep[j]>0)||(cutinfo->nY[j]>1))
			//	continue;

			//if ((j0==j)||(j0==j1)||(cutinfo->FDep[j0]>0)||(cutinfo->nY[j0]>1))
			//	continue;

			//if (fabs(cutinfo->p[j1]-(cutinfo->p[j]+cutinfo->p[j0]))>0.001)
			//	continue;

			//if (cutinfo->u[j1] < (cutinfo->u[j]+cutinfo->u[j0]) - 0.001)
			//	continue;

			////if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
			////	continue;

			for (i1 = 0; i1 < m; i1++)
			{
				k1 = i1 * n + j1;
				if (x[k1] < 0.001)
					continue;

				for (i = 0; i < i1 - 1; i++)
				{
					k = i * n + j;
					k0 = i * n + j0;
					cutvio = x[k] + x[k0];

					if ((cutvio + x[k1]) > 2.001)
					{
						//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
						//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
						//for (i1=0; i1<KP.n; i1++)
						//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

						kk = 0;
						rhs[0] = 2.0;
						beg[0] = kk;

						ind[kk] = k1;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k0;
						val[kk] = 1.0;
						kk++;

						beg[1] = kk;
						cutinfo->nz = kk;

						return 1;
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_2ItemRev(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, h, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag2IRev[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (h = 0; h < cutinfo->Flag2IRev[j1]; h++)
		{
			j = cutinfo->Tupla2IRev[j1][h].j1;
			j0 = cutinfo->Tupla2IRev[j1][h].j2;

			//if ((j==j1)||(cutinfo->FDep[j]>0)||(cutinfo->nY[j]>1))
			//	continue;

			//if ((j0==j)||(j0==j1)||(cutinfo->FDep[j0]>0)||(cutinfo->nY[j0]>1))
			//	continue;

			//if (fabs(cutinfo->p[j1]-(cutinfo->p[j]+cutinfo->p[j0]))>0.001)
			//	continue;

			//if (cutinfo->u[j1] > (cutinfo->u[j]+cutinfo->u[j0]) - 0.001)
			//	continue;

			////if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
			////	continue;

			for (i1 = 0; i1 < m - 1; i1++)
			{
				k1 = i1 * n + j1;
				if (x[k1] < 0.001)
					continue;

				for (i = i1 + 1; i < m; i++)
				{
					k = i * n + j;
					k0 = i * n + j0;
					cutvio = x[k] + x[k0];

					if ((cutvio + x[k1]) > 2.001)
					{
						//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
						//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
						//for (i1=0; i1<KP.n; i1++)
						//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

						kk = 0;
						rhs[0] = 2.0;
						beg[0] = kk;

						ind[kk] = k1;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k0;
						val[kk] = 1.0;
						kk++;

						beg[1] = kk;
						cutinfo->nz = kk;

						return 1;
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_3Item(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, h, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0, j2, k2;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag3I[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (h = 0; h < cutinfo->Flag3I[j1]; h++)
		{
			j = cutinfo->Tupla3I[j1][h].j1;
			j0 = cutinfo->Tupla3I[j1][h].j2;
			j2 = cutinfo->Tupla3I[j1][h].j3;

			//if ((j==j1)||(cutinfo->FDep[j]>0)||(cutinfo->nY[j]>1))
			//	continue;

			//if ((j0==j)||(j0==j1)||(cutinfo->FDep[j0]>0)||(cutinfo->nY[j0]>1))
			//	continue;

			//if ((j2==j)||(j2==j1)||(j2==j0)||(cutinfo->FDep[j2]>0)||(cutinfo->nY[j2]>1))
			//	continue;

			//if (fabs(cutinfo->p[j1]-(cutinfo->p[j]+cutinfo->p[j0]+cutinfo->p[j2]))>0.001)
			//	continue;

			//if (cutinfo->u[j1] < (cutinfo->u[j]+cutinfo->u[j0]+cutinfo->u[j2]) - 0.001)
			//	continue;

			////if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
			////	continue;

			for (i1 = 0; i1 < m; i1++)
			{
				k1 = i1 * n + j1;
				if (x[k1] < 0.001)
					continue;

				for (i = 0; i < i1 - 1; i++)
				{
					k = i * n + j;
					k0 = i * n + j0;
					k2 = i * n + j2;
					cutvio = x[k] + x[k0] + x[k2];

					if ((cutvio + x[k1]) > 3.001)
					{
						//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
						//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
						//for (i1=0; i1<KP.n; i1++)
						//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

						kk = 0;
						rhs[0] = 3.0;
						beg[0] = kk;

						ind[kk] = k1;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k0;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k2;
						val[kk] = 1.0;
						kk++;

						beg[1] = kk;
						cutinfo->nz = kk;

						return 1;
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_3ItemRev(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, h, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0, j2, k2;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag3IRev[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (h = 0; h < cutinfo->Flag3IRev[j1]; h++)
		{
			j = cutinfo->Tupla3IRev[j1][h].j1;
			j0 = cutinfo->Tupla3IRev[j1][h].j2;
			j2 = cutinfo->Tupla3IRev[j1][h].j3;

			//if ((j==j1)||(cutinfo->FDep[j]>0)||(cutinfo->nY[j]>1))
			//	continue;

			//if ((j0==j)||(j0==j1)||(cutinfo->FDep[j0]>0)||(cutinfo->nY[j0]>1))
			//	continue;

			//if ((j2==j)||(j2==j1)||(j2==j0)||(cutinfo->FDep[j2]>0)||(cutinfo->nY[j2]>1))
			//	continue;

			//if (fabs(cutinfo->p[j1]-(cutinfo->p[j]+cutinfo->p[j0]+cutinfo->p[j2]))>0.001)
			//	continue;

			//if (cutinfo->u[j1] > (cutinfo->u[j]+cutinfo->u[j0]+cutinfo->u[j2]) - 0.001)
			//	continue;

			////if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
			////	continue;

			for (i1 = 0; i1 < m - 1; i1++)
			{
				k1 = i1 * n + j1;
				if (x[k1] < 0.001)
					continue;

				for (i = i1 + 1; i < m; i++)
				{
					k = i * n + j;
					k0 = i * n + j0;
					k2 = i * n + j2;
					cutvio = x[k] + x[k0] + x[k2];

					if ((cutvio + x[k1]) > 3.001)
					{
						//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
						//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
						//for (i1=0; i1<KP.n; i1++)
						//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

						kk = 0;
						rhs[0] = 3.0;
						beg[0] = kk;

						ind[kk] = k1;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k0;
						val[kk] = 1.0;
						kk++;

						ind[kk] = k2;
						val[kk] = 1.0;
						kk++;

						beg[1] = kk;
						cutinfo->nz = kk;

						return 1;
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Cover Inequalities
//-----------------------------------------------------------------------------
//
int Cover_Inequalities(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;
	long     konta, cover;
	long     ww, minw, maxw;
	double   lhs, maxlhs;
	//Knapsack KP;

	// Setup Knapsack data structure
	//KP.Malloc(n);

	// Check a cover for each knapsack constraint
	for (i = 0; i < m; i++)
	{
		cutinfo->KP.n = 0;
		cutinfo->KP.nf = 0;
		cutinfo->KP.zf = 0;
		//cutinfo->KP.W = int(cutinfo->pmax[i]+0.001);
		cutinfo->KP.W = int(cutinfo->pmax[i] * 10 + 0.001);

		for (j = 0; j < n; j++)
		{
			cutinfo->KP.fcover[j] = 0;
			k = i * n + j;
			if (x[k] < 0.001) continue;

			if (x[k] > 0.99999999)
			{
				cutinfo->KP.fcover[j] = 1;
				cutinfo->KP.wf[cutinfo->KP.nf] = int(cutinfo->p[j] * 10 + 0.001);
				cutinfo->KP.indf[cutinfo->KP.nf] = j;
				cutinfo->KP.zf = cutinfo->KP.zf + 1;
				cutinfo->KP.W = cutinfo->KP.W - cutinfo->KP.wf[cutinfo->KP.nf];
				cutinfo->KP.nf++;
			}
			else
			{
				//cutinfo->KP.w[cutinfo->KP.n] = int(cutinfo->p[j]+0.001);
				cutinfo->KP.w[cutinfo->KP.n] = int(cutinfo->p[j] * 10 + 0.001);
				cutinfo->KP.p[cutinfo->KP.n] = x[k];
				cutinfo->KP.ind[cutinfo->KP.n] = j;
				cutinfo->KP.n++;
			}
		}

		if (cutinfo->KP.n > 0)
		{
#ifdef ReMap
			cutinfo->KP.SolveReMap();
#else
			cutinfo->KP.Solve();
#endif
			cutinfo->KP.zopt = cutinfo->KP.zopt + cutinfo->KP.zf;
		}
		else
			cutinfo->KP.zopt = cutinfo->KP.zf;

		konta = cutinfo->KP.nf;
		for (i1 = 0; i1 < cutinfo->KP.n; i1++)
		{
			if (cutinfo->KP.x[i1] > 0)
			{
				cutinfo->KP.fcover[cutinfo->KP.ind[i1]] = 1;
				konta++;
			}
		}

		// Check if a violated minimal cover exists
		cover = -1;
		//minw=100000000;
		maxlhs = 0.0;
		for (i1 = 0; i1 < cutinfo->KP.n; i1++)
		{
			if (cutinfo->KP.x[i1] > 0) continue;
			if (((double)cutinfo->KP.zopt + cutinfo->KP.p[i1]) > ((double)konta + 0.0001))
			{
				lhs = ((double)cutinfo->KP.zopt + cutinfo->KP.p[i1]);
				if (lhs > maxlhs)
				{
					lhs = maxlhs;
					cover = cutinfo->KP.ind[i1];
				}
				//if (minw>cutinfo->KP.w[i1])
				//{
				//	minw=cutinfo->KP.w[i1];
				//	cover=cutinfo->KP.ind[i1];
				//}
			}
		}

		if (cover > -1)
		{
			cutinfo->KP.fcover[cover] = 1;

			//printf("\n Knapsack  W=%d  zopt=%lf  zf=%lf\n",cutinfo->KP.W,cutinfo->KP.zopt,cutinfo->KP.zf);
			////for (i1=0; i1<cutinfo->KP.n; i1++)
			////	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,cutinfo->KP.ind[i1],cutinfo->KP.p[i1],KP.w[i1],cutinfo->KP.x[i1]);

			// Try to extend the cover
			maxw = 0;
			for (j = 0; j < n; j++)
			{
				if (cutinfo->KP.fcover[j])
				{
					ww = int(cutinfo->p[j] * 10 + 0.001);
					if (ww > maxw)
						maxw = ww;
				}
			}

			for (j = 0; j < n; j++)
			{
				cutinfo->KP.wf[j] = 1;
				if (cutinfo->KP.fcover[j]) continue;

				ww = int(cutinfo->p[j] * 10 + 0.001);
				if (ww > maxw)
					cutinfo->KP.fcover[j] = 1;
			}

			// Try to lift the cover inequality found
			if (cutinfo->Lifting)
			{
				cutinfo->KP.nf = 0;
				for (j = 0; j < n; j++)
				{
					if (cutinfo->KP.fcover[j] > 0)
					{
						cutinfo->KP.wf[j] = 1;
						cutinfo->KP.indf[cutinfo->KP.nf] = j;
						cutinfo->KP.nf++;
					}
				}
				for (i1 = 0; i1 < cutinfo->KP.nf; i1++)
				{
					j1 = cutinfo->KP.indf[i1];
					cutinfo->KP.n = 0;
					cutinfo->KP.W = int(cutinfo->pmax[i] * 10 - cutinfo->p[j1] * 10 + 0.001);
					for (j = 0; j < n; j++)
					{
						if ((cutinfo->KP.fcover[j] > 0) && (j != j1))
						{
							cutinfo->KP.w[cutinfo->KP.n] = int(cutinfo->p[j] * 10 + 0.001);
							cutinfo->KP.p[cutinfo->KP.n] = (double)cutinfo->KP.wf[j];
							cutinfo->KP.ind[cutinfo->KP.n] = j;
							cutinfo->KP.n++;
						}
					}
#ifdef ReMap
					cutinfo->KP.SolveReMap();
#else
					cutinfo->KP.Solve();
#endif
					if (((double)konta - cutinfo->KP.zopt) > 1.001)
					{
						cutinfo->KP.wf[j1] = int((double)konta - cutinfo->KP.zopt + 0.001);
					}
				}
			}

			// Save the cover inequality
			k1 = 0;
			cutinfo->rhs[0] = (double)konta;
			cutinfo->beg[0] = k1;
			for (j = 0; j < n; j++)
			{
				if (cutinfo->KP.fcover[j])
				{
					k = i * n + j;
					//ww=int(cutinfo->p[j]*10+0.001);
					//printf("%d)  x(%d)=p=%lf   w=%d    x=1 *\n",k1,j,x[k],ww);
					cutinfo->ind[k1] = k;
					//cutinfo->val[k1]=1.0; // Prima 1.0
					cutinfo->val[k1] = (double)cutinfo->KP.wf[j];
					k1++;
				}
			}
			//for (i1=0; i1<cutinfo->KP.nf; i1++)
			//{
			//	k = i*n+cutinfo->KP.indf[i1];
			//	//printf("%d)  x(%d)=p=%lf   w=%d    x=1 *\n",k1,cutinfo->KP.indf[i1],x[k],cutinfo->KP.wf[i1]);
			//	cutinfo->ind[k1]=k;
			//	cutinfo->val[k1]=1.0;
			//	k1++;
			//}
			//for (i1=0; i1<cutinfo->KP.n; i1++)
			//	if ((cover==i1)||(cutinfo->KP.x[i1]>0))
			//	{
			//		k = i*n+cutinfo->KP.ind[i1];
			//		//printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",k1,cutinfo->KP.ind[i1],x[k],cutinfo->KP.w[i1],cutinfo->KP.x[i1]);
			//		cutinfo->ind[k1]=k;
			//		cutinfo->val[k1]=1.0;
			//		k1++;
			//	}
			cutinfo->beg[1] = k1;
			cutinfo->nz = k1;

			return 1;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Cover Inequalities using Cplex
//-----------------------------------------------------------------------------
//
int Cover_InequalitiesCplex(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;

	long    konta, cover;
	double  ww, minw, maxw;
	double  lhs, maxlhs;
	double  ub;
	double  zopt;
	long* fcover;
	double* ccover;

	fcover = new long[n];
	ccover = new double[n];
	for (j = 0; j < n; j++)
	{
		ind[j] = j;
	}

	// Check a cover for each knapsack constraint
	for (i = 0; i < m; i++)
	{
		// RHS: Knapsack Capacity
		status = CPXchgrhs(cutinfo->env, cutinfo->lpkp, 1, ind, &(cutinfo->pmax[i]));
		status = CPXchgobj(cutinfo->env, cutinfo->lpkp, n, ind, x + i * n);
		status = CPXmipopt(cutinfo->env, cutinfo->lpkp);
		status = CPXgetobjval(cutinfo->env, cutinfo->lpkp, &zopt);
		status = CPXgetx(cutinfo->env, cutinfo->lpkp, val, 0, n - 1);

		konta = 0;
		for (j = 0; j < n; j++)
		{
			if (val[j] > 0)
			{
				fcover[j] = 1;
				ccover[j] = 1.0;
				konta++;
			}
			else
			{
				fcover[j] = 0;
				ccover[j] = 0.0;
			}
		}

		// Check if a violated minimal cover exists
		cover = -1;
		maxlhs = 0.0;
		for (j = 0; j < n; j++)
		{
			if (val[j] > 0) continue;
			if ((zopt + x[i * n + j]) > ((double)konta + 0.0001))
			{
				lhs = zopt + x[i * n + j];
				if (lhs > maxlhs)
				{
					lhs = maxlhs;
					cover = j;
				}
			}
		}

		if (cover > -1)
		{
			fcover[cover] = 1;
			ccover[cover] = 1.0;

			//printf("\n Knapsack  W=%d  zopt=%lf  zf=%lf\n",KP.W,KP.zopt,KP.zf);
			////for (i1=0; i1<KP.n; i1++)
			////	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

			// Try to extend the cover
			maxw = 0.0;
			for (j = 0; j < n; j++)
			{
				if (fcover[j])
				{
					if (cutinfo->p[j] > maxw)
						maxw = cutinfo->p[j];
				}
			}

			for (j = 0; j < n; j++)
			{
				if (fcover[j]) continue;

				if (cutinfo->p[j] > maxw)
				{
					fcover[j] = 1;
					ccover[j] = 1.0;
				}
			}

			// Try to lift the cover inequality found
			if (cutinfo->Lifting)
			{
				for (j = 0; j < n; j++)
				{
					if (fcover[j] > 0)
					{
						ub = 1.0;
						status = CPXchgbds(cutinfo->env, cutinfo->lpkp, 1, &j, "L", &ub);
						status = CPXmipopt(cutinfo->env, cutinfo->lpkp);
						status = CPXgetobjval(cutinfo->env, cutinfo->lpkp, &zopt);

						if (((double)konta - zopt) > 1.001)
							ccover[j] = int((double)konta - zopt + 0.001);

						ub = 0.0;
						status = CPXchgbds(cutinfo->env, cutinfo->lpkp, 1, &j, "L", &ub);
					}
				}
			}

			// Save the cover inequality
			k1 = 0;
			cutinfo->rhs[0] = (double)konta;
			cutinfo->beg[0] = k1;
			for (j = 0; j < n; j++)
			{
				if (fcover[j])
				{
					k = i * n + j;
					cutinfo->ind[k1] = k;
					cutinfo->val[k1] = (double)ccover[j];
					k1++;
				}
			}

			cutinfo->beg[1] = k1;
			cutinfo->nz = k1;

			return 1;
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Cover Inequalities
//-----------------------------------------------------------------------------
//
int Cover_InequalitiesOld(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1;
	long     konta, cover;
	long     ww, minw, maxw;
	double   lhs, maxlhs;
	Knapsack KP;

	maxw = 0;
	for (i = 0; i < m; i++)
	{
		if (maxw<int(cutinfo->pmax[i] * 10 + 0.001))
			maxw = int(cutinfo->pmax[i] * 10 + 0.001);
	}

	// Setup Knapsack data structure
	KP.Malloc(n, maxw);

	// Check a cover for each knapsack constraint
	for (i = 0; i < m; i++)
	{
		KP.n = 0;
		KP.nf = 0;
		KP.zf = 0;
		//KP.W = int(cutinfo->pmax[i]+0.001);
		KP.W = int(cutinfo->pmax[i] * 10 + 0.001);

		for (j = 0; j < n; j++)
		{
			KP.fcover[j] = 0;
			k = i * n + j;
			if (x[k] < 0.001) continue;

			if (x[k] > 0.99999999)
			{
				KP.fcover[j] = 1;
				KP.wf[KP.nf] = int(cutinfo->p[j] * 10 + 0.001);
				KP.indf[KP.nf] = j;
				KP.zf = KP.zf + 1;
				KP.W = KP.W - KP.wf[KP.nf];
				KP.nf++;
			}
			else
			{
				//KP.w[KP.n] = int(cutinfo->p[j]+0.001);
				KP.w[KP.n] = int(cutinfo->p[j] * 10 + 0.001);
				KP.p[KP.n] = x[k];
				KP.ind[KP.n] = j;
				KP.n++;
			}
		}

		if (KP.n > 0)
		{
#ifdef ReMap
			KP.SolveReMap();
#else
			KP.Solve();
#endif
			KP.zopt = KP.zopt + KP.zf;
		}
		else
			KP.zopt = KP.zf;

		konta = KP.nf;
		for (i1 = 0; i1 < KP.n; i1++)
		{
			if (KP.x[i1] > 0)
			{
				KP.fcover[KP.ind[i1]] = 1;
				konta++;
			}
		}

		// Check if a violated minimal cover exists
		cover = -1;
		//minw=100000000;
		maxlhs = 0.0;
		for (i1 = 0; i1 < KP.n; i1++)
		{
			if (KP.x[i1] > 0) continue;
			if (((double)KP.zopt + KP.p[i1]) > ((double)konta + 0.0001))
			{
				lhs = ((double)KP.zopt + KP.p[i1]);
				if (lhs > maxlhs)
				{
					lhs = maxlhs;
					cover = KP.ind[i1];
				}
				//if (minw>KP.w[i1])
				//{
				//	minw=KP.w[i1];
				//	cover=KP.ind[i1];
				//}
			}
		}

		if (cover > -1)
		{
			KP.fcover[cover] = 1;

			//printf("\n Knapsack  W=%d  zopt=%lf  zf=%lf\n",KP.W,KP.zopt,KP.zf);
			////for (i1=0; i1<KP.n; i1++)
			////	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

			// Try to extend the cover
			maxw = 0;
			for (j = 0; j < n; j++)
			{
				if (KP.fcover[j])
				{
					ww = int(cutinfo->p[j] * 10 + 0.001);
					if (ww > maxw)
						maxw = ww;
				}
			}

			for (j = 0; j < n; j++)
			{
				KP.wf[j] = 1;
				if (KP.fcover[j]) continue;

				ww = int(cutinfo->p[j] * 10 + 0.001);
				if (ww > maxw)
					KP.fcover[j] = 1;
			}

			// Try to lift the cover inequality found
			if (cutinfo->Lifting)
			{
				KP.nf = 0;
				for (j = 0; j < n; j++)
				{
					if (KP.fcover[j] > 0)
					{
						KP.wf[j] = 1;
						KP.indf[KP.nf] = j;
						KP.nf++;
					}
				}
				for (i1 = 0; i1 < KP.nf; i1++)
				{
					j1 = KP.indf[i1];
					KP.n = 0;
					KP.W = int(cutinfo->pmax[i] * 10 - cutinfo->p[j1] * 10 + 0.001);
					for (j = 0; j < n; j++)
					{
						if ((KP.fcover[j] > 0) && (j != j1))
						{
							KP.w[KP.n] = int(cutinfo->p[j] * 10 + 0.001);
							KP.p[KP.n] = (double)KP.wf[j];
							KP.ind[KP.n] = j;
							KP.n++;
						}
					}
#ifdef ReMap
					KP.SolveReMap();
#else
					KP.Solve();
#endif
					if (((double)konta - KP.zopt) > 1.001)
					{
						KP.wf[j1] = int((double)konta - KP.zopt + 0.001);
					}
				}
			}

			// Save the cover inequality
			k1 = 0;
			cutinfo->rhs[0] = (double)konta;
			cutinfo->beg[0] = k1;
			for (j = 0; j < n; j++)
			{
				if (KP.fcover[j])
				{
					k = i * n + j;
					//ww=int(cutinfo->p[j]*10+0.001);
					//printf("%d)  x(%d)=p=%lf   w=%d    x=1 *\n",k1,j,x[k],ww);
					cutinfo->ind[k1] = k;
					//cutinfo->val[k1]=1.0; // Prima 1.0
					cutinfo->val[k1] = (double)KP.wf[j];
					k1++;
				}
			}
			//for (i1=0; i1<KP.nf; i1++)
			//{
			//	k = i*n+KP.indf[i1];
			//	//printf("%d)  x(%d)=p=%lf   w=%d    x=1 *\n",k1,KP.indf[i1],x[k],KP.wf[i1]);
			//	cutinfo->ind[k1]=k;
			//	cutinfo->val[k1]=1.0;
			//	k1++;
			//}
			//for (i1=0; i1<KP.n; i1++)
			//	if ((cover==i1)||(KP.x[i1]>0))
			//	{
			//		k = i*n+KP.ind[i1];
			//		//printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",k1,KP.ind[i1],x[k],KP.w[i1],KP.x[i1]);
			//		cutinfo->ind[k1]=k;
			//		cutinfo->val[k1]=1.0;
			//		k1++;
			//	}
			cutinfo->beg[1] = k1;
			cutinfo->nz = k1;

			return 1;
		}
	}

	return 0;
}


int CalcolaFlag2I(CUTINFOptr cutinfo)
{
	long j, j1, j0;

	for (j1 = 0; j1 < cutinfo->n; j1++)
	{
		cutinfo->Flag2I[j1] = 0;
		cutinfo->Tupla2I[j1] = NULL;

		if ((cutinfo->FDepA[j1] > 0) || (cutinfo->nY[j1] > 1))
			//if ((cutinfo->nY[j1]>1))
			continue;

		for (j = 0; j < cutinfo->n - 1; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0])) > 0.001)
					continue;

				if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0]) - 0.001)
					continue;

				cutinfo->Flag2I[j1]++;

			}
		}

		if (cutinfo->Flag2I[j1] == 0)
			continue;

		cutinfo->Tupla2I[j1] = new struct tupla[cutinfo->Flag2I[j1]];

		cutinfo->Flag2I[j1] = 0;
		for (j = 0; j < cutinfo->n - 1; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0])) > 0.001)
					continue;

				if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0]) - 0.001)
					continue;

				cutinfo->Tupla2I[j1][cutinfo->Flag2I[j1]].j1 = j;
				cutinfo->Tupla2I[j1][cutinfo->Flag2I[j1]].j2 = j0;

				cutinfo->Flag2I[j1]++;

			}
		}
	}

	return 0;
}


int CalcolaFlag3I(CUTINFOptr cutinfo)
{
	long j, j1, j0, j2;

	for (j1 = 0; j1 < cutinfo->n; j1++)
	{
		cutinfo->Flag3I[j1] = 0;
		cutinfo->Tupla3I[j1] = NULL;

		if ((cutinfo->FDepA[j1] > 0) || (cutinfo->nY[j1] > 1))
			//if ((cutinfo->nY[j1]>1))
			continue;

		for (j = 0; j < cutinfo->n - 2; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n - 1; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				for (j2 = j0 + 1; j2 < cutinfo->n; j2++)
				{
					if ((j2 == j) || (j2 == j1) || (j2 == j0) || (cutinfo->FDepP[j2] > 0) || (cutinfo->nY[j2] > 1))
						continue;

					if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0] + cutinfo->p[j2])) > 0.001)
						continue;

					if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0] + cutinfo->u[j2]) - 0.001)
						continue;

					cutinfo->Flag3I[j1]++;

				}
			}
		}

		if (cutinfo->Flag3I[j1] == 0)
			continue;

		cutinfo->Tupla3I[j1] = new struct tupla[cutinfo->Flag3I[j1]];

		cutinfo->Flag3I[j1] = 0;
		for (j = 0; j < cutinfo->n - 2; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n - 1; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				for (j2 = j0 + 1; j2 < cutinfo->n; j2++)
				{
					if ((j2 == j) || (j2 == j1) || (j2 == j0) || (cutinfo->FDepP[j2] > 0) || (cutinfo->nY[j2] > 1))
						continue;

					if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0] + cutinfo->p[j2])) > 0.001)
						continue;

					if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0] + cutinfo->u[j2]) - 0.001)
						continue;

					cutinfo->Tupla3I[j1][cutinfo->Flag3I[j1]].j1 = j;
					cutinfo->Tupla3I[j1][cutinfo->Flag3I[j1]].j2 = j0;
					cutinfo->Tupla3I[j1][cutinfo->Flag3I[j1]].j3 = j2;

					cutinfo->Flag3I[j1]++;

				}
			}
		}
	}

	return 0;
}


int CalcolaFlag2IRev(CUTINFOptr cutinfo)
{
	long j, j1, j0;

	for (j1 = 0; j1 < cutinfo->n; j1++)
	{
		cutinfo->Flag2IRev[j1] = 0;
		cutinfo->Tupla2IRev[j1] = NULL;

		if ((cutinfo->FDepP[j1] > 0) || (cutinfo->nY[j1] > 1))
			//if ((cutinfo->nY[j1]>1))
			continue;

		for (j = 0; j < cutinfo->n - 1; j++)
		{
			if ((j == j1) || (cutinfo->FDepA[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepA[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0])) > 0.001)
					continue;

				if (cutinfo->u[j1] > (cutinfo->u[j] + cutinfo->u[j0]) - 0.001)
					continue;

				cutinfo->Flag2IRev[j1]++;

			}
		}

		if (cutinfo->Flag2IRev[j1] == 0)
			continue;

		cutinfo->Tupla2IRev[j1] = new struct tupla[cutinfo->Flag2IRev[j1]];

		cutinfo->Flag2IRev[j1] = 0;
		for (j = 0; j < cutinfo->n - 1; j++)
		{
			if ((j == j1) || (cutinfo->FDepA[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepA[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0])) > 0.001)
					continue;

				if (cutinfo->u[j1] > (cutinfo->u[j] + cutinfo->u[j0]) - 0.001)
					continue;

				cutinfo->Tupla2IRev[j1][cutinfo->Flag2IRev[j1]].j1 = j;
				cutinfo->Tupla2IRev[j1][cutinfo->Flag2IRev[j1]].j2 = j0;

				cutinfo->Flag2IRev[j1]++;

			}
		}
	}

	return 0;
}


int CalcolaFlag3IRev(CUTINFOptr cutinfo)
{
	long j, j1, j0, j2;

	for (j1 = 0; j1 < cutinfo->n; j1++)
	{
		cutinfo->Flag3IRev[j1] = 0;
		cutinfo->Tupla3IRev[j1] = NULL;

		if ((cutinfo->FDepP[j1] > 0) || (cutinfo->nY[j1] > 1))
			//if ((cutinfo->nY[j1]>1))
			continue;

		for (j = 0; j < cutinfo->n - 2; j++)
		{
			if ((j == j1) || (cutinfo->FDepA[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n - 1; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepA[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				for (j2 = j0 + 1; j2 < cutinfo->n; j2++)
				{
					if ((j2 == j) || (j2 == j1) || (j2 == j0) || (cutinfo->FDepA[j2] > 0) || (cutinfo->nY[j2] > 1))
						continue;

					if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0] + cutinfo->p[j2])) > 0.001)
						continue;

					if (cutinfo->u[j1] > (cutinfo->u[j] + cutinfo->u[j0] + cutinfo->u[j2]) - 0.001)
						continue;

					cutinfo->Flag3IRev[j1]++;

				}
			}
		}

		if (cutinfo->Flag3IRev[j1] == 0)
			continue;

		cutinfo->Tupla3IRev[j1] = new struct tupla[cutinfo->Flag3IRev[j1]];

		cutinfo->Flag3IRev[j1] = 0;
		for (j = 0; j < cutinfo->n - 2; j++)
		{
			if ((j == j1) || (cutinfo->FDepA[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < cutinfo->n - 1; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepA[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				for (j2 = j0 + 1; j2 < cutinfo->n; j2++)
				{
					if ((j2 == j) || (j2 == j1) || (j2 == j0) || (cutinfo->FDepA[j2] > 0) || (cutinfo->nY[j2] > 1))
						continue;

					if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0] + cutinfo->p[j2])) > 0.001)
						continue;

					if (cutinfo->u[j1] > (cutinfo->u[j] + cutinfo->u[j0] + cutinfo->u[j2]) - 0.001)
						continue;

					cutinfo->Tupla3IRev[j1][cutinfo->Flag3IRev[j1]].j1 = j;
					cutinfo->Tupla3IRev[j1][cutinfo->Flag3IRev[j1]].j2 = j0;
					cutinfo->Tupla3IRev[j1][cutinfo->Flag3IRev[j1]].j3 = j2;

					cutinfo->Flag3IRev[j1]++;

				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_2Item_Old(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag2I[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (j = 0; j < n - 1; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < n; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0])) > 0.001)
					continue;

				if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0]) - 0.001)
					continue;

				//if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
				//	continue;

				for (i1 = 0; i1 < m; i1++)
				{
					k1 = i1 * n + j1;
					if (x[k1] < 0.001)
						continue;

					for (i = 0; i < i1 - 1; i++)
					{
						k = i * n + j;
						k0 = i * n + j0;
						cutvio = x[k] + x[k0];

						if ((cutvio + x[k1]) > 2.001)
						{
							//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
							//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
							//for (i1=0; i1<KP.n; i1++)
							//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

							kk = 0;
							rhs[0] = 2.0;
							beg[0] = kk;

							ind[kk] = k1;
							val[kk] = 1.0;
							kk++;

							ind[kk] = k;
							val[kk] = 1.0;
							kk++;

							ind[kk] = k0;
							val[kk] = 1.0;
							kk++;

							beg[1] = kk;
							cutinfo->nz = kk;

							return 1;
						}
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Cutting plane procedure for Precedence constraints
//-----------------------------------------------------------------------------
//
int Precedence_Inequalities_3Item_Old(CUTINFOptr cutinfo)
{
	int status = 0;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double* x = cutinfo->x;
	int* beg = cutinfo->beg;
	int* ind = cutinfo->ind;
	double* val = cutinfo->val;
	double* rhs = cutinfo->rhs;
	long     n = cutinfo->n;
	long     m = cutinfo->m;
	double   cutvio;
	int      i, j, k, cutnz;
	long     i1, j1, k1, kk;
	long     j0, k0, j2, k2;

	// Precedence search
	for (j1 = 0; j1 < n; j1++)
	{
		//if ((cutinfo->FDep[j1]>0)||(cutinfo->nY[j1]>1))
		if ((cutinfo->Flag3I[j1] == 0) || (cutinfo->nY[j1] > 1))
			continue;

		for (j = 0; j < n - 2; j++)
		{
			if ((j == j1) || (cutinfo->FDepP[j] > 0) || (cutinfo->nY[j] > 1))
				continue;

			for (j0 = j + 1; j0 < n - 1; j0++)
			{
				if ((j0 == j) || (j0 == j1) || (cutinfo->FDepP[j0] > 0) || (cutinfo->nY[j0] > 1))
					continue;

				for (j2 = j0 + 1; j2 < n; j2++)
				{
					if ((j2 == j) || (j2 == j1) || (j2 == j0) || (cutinfo->FDepP[j2] > 0) || (cutinfo->nY[j2] > 1))
						continue;

					if (fabs(cutinfo->p[j1] - (cutinfo->p[j] + cutinfo->p[j0] + cutinfo->p[j2])) > 0.001)
						continue;

					if (cutinfo->u[j1] < (cutinfo->u[j] + cutinfo->u[j0] + cutinfo->u[j2]) - 0.001)
						continue;

					//if ((fabs(cutinfo->u[j1]-(cutinfo->u[j]+cutinfo->u[j0]))<0.001)&&(j1>j)&&(j1>j0))
					//	continue;

					for (i1 = 0; i1 < m; i1++)
					{
						k1 = i1 * n + j1;
						if (x[k1] < 0.001)
							continue;

						for (i = 0; i < i1 - 1; i++)
						{
							k = i * n + j;
							k0 = i * n + j0;
							k2 = i * n + j2;
							cutvio = x[k] + x[k0] + x[k2];

							if ((cutvio + x[k1]) > 3.001)
							{
								//printf("\n Precedence:  (j1,u,p)=(%3d,%5.1lf,%5.1lf)  (j,u,p)=(%3d,%5.1lf,%5.1lf)  Sprint: %d\n",
								//	   j1,cutinfo->u[j1],cutinfo->p[j1],j,cutinfo->u[j],cutinfo->p[j],i1);
								//for (i1=0; i1<KP.n; i1++)
								//	printf("%d)  x(%d)=p=%lf   w=%d    x=%d \n",i1,KP.ind[i1],KP.p[i1],KP.w[i1],KP.x[i1]);

								kk = 0;
								rhs[0] = 3.0;
								beg[0] = kk;

								ind[kk] = k1;
								val[kk] = 1.0;
								kk++;

								ind[kk] = k;
								val[kk] = 1.0;
								kk++;

								ind[kk] = k0;
								val[kk] = 1.0;
								kk++;

								ind[kk] = k2;
								val[kk] = 1.0;
								kk++;

								beg[1] = kk;
								cutinfo->nz = kk;

								return 1;
							}
						}
					}
				}
			}
		}
	}

	return 0;
}


//-----------------------------------------------------------------------------
// Setup the knapsack subproblem
//-----------------------------------------------------------------------------
//
int  CplexObj::SetupKPDP(CUTINFOptr cutinfo)
{
	long i, maxw;

	maxw = 0;
	for (i = 0; i < cutinfo->m; i++)
		if (maxw<int(cutinfo->pmax[i] * 10 + 0.001))
			maxw = int(cutinfo->pmax[i] * 10 + 0.001);

	cutinfo->KP.Malloc(cutinfo->n, maxw);

	return 0;
}

//-----------------------------------------------------------------------------
// Setup the knapsack subproblem
//-----------------------------------------------------------------------------
//
int  CplexObj::SetupKPMIP(CUTINFOptr cutinfo)
{
	int status = 0;
	int ncols;
	int nrows;
	int objsen;
	double objval;
	double* obj;
	double* rhs;
	char* sense;
	int* matbeg;
	int* matcnt;
	int* matind;
	double* matval;
	double* lb;
	double* ub;
	double* xt;
	int* indices;
	char* xctype;
	double* rngval = NULL;

	// Auxiliary variables
	long kc, knz;
	long i;

	if (env == NULL)
	{
		printf("\nError... env not alocated...\n");
		return 2;
	}

	cutinfo->lpkp = CPXcreateprob(env, &status, "lpkp");
	if (status)
	{
		printf("\nError... CPXcreateprob...\n");
		lp = NULL;
		return 1;
	}

	// Variables
	ncols = cutinfo->n;

	// Constraints
	nrows = 1;

	// Maximization
	objsen = -1;

	obj = new double[ncols];
	rhs = new double[nrows];
	sense = new char[nrows];
	matbeg = new int[ncols + 4];
	matcnt = new int[ncols + 4];
	matind = new int[ncols];
	matval = new double[ncols];
	lb = new double[ncols];
	ub = new double[ncols];
	xt = new double[ncols];
	indices = new int[ncols];
	xctype = new char[ncols];

	kc = 0;
	knz = 0;
	for (i = 0;i < ncols;i++)
	{
		// Objective function
		obj[kc] = 1.0;

		// Column
		matbeg[kc] = knz;
		matcnt[kc] = 1;

		matind[knz] = 0;
		matval[knz] = cutinfo->p[i];
		knz++;

		// Lower and upper bound
		lb[kc] = 0.;
		ub[kc] = 1.;
		indices[kc] = kc;
		xctype[kc] = 'B';

		kc++;
	}

	matbeg[kc] = knz;

	// Right Hand Sides
	rhs[0] = cutinfo->pmax[0];
	sense[0] = 'L';

	rngval = NULL;

	// Load the data into Cplex
	status = CPXcopylp(env, cutinfo->lpkp, ncols, nrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, rngval);

	status = CPXchgprobtype(env, cutinfo->lpkp, CPXPROB_MILP);
	status = CPXchgctype(env, cutinfo->lpkp, ncols, indices, xctype);
	//status = CPXwriteprob(env,cutinfo->lpkp,"mysep.lp","LP");

	// Free memory
	delete[] obj;
	delete[] rhs;
	delete[] sense;
	delete[] matbeg;
	delete[] matcnt;
	delete[] matind;
	delete[] matval;
	delete[] lb;
	delete[] ub;
	delete[] xt;
	delete[] indices;
	delete[] xctype;

	return 0;
}
