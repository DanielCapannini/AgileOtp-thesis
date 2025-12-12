//-----------------------------------------------------------------------------
//  File: CplexObj.h                                                       
//                                                                           
//  Project: Cplex solver interface (header)                   
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 24.08.2011                                                  
//-----------------------------------------------------------------------------
//
#include "cplex.h"
#include "knapsack.h"

struct tupla {
	long j1,j2,j3;
};

struct cutinfo {
	int FlagLP;   // FlagLP=0 if we are solving the original model; FlagLP=1 if we are solving the knapsack subproblem; 
	CPXENVptr env;  // Pointer to Cplex Environment 
	CPXLPptr lp;
	CPXLPptr lpkp;
	int      numcols;
	int      numtoadd;
	int      maxcuts;
	int      num;
	int      nz;
	double   *x;
	int      *beg;
	int      *ind; 
	double   *val;
	double   *rhs;
	int      *Flag2I;
	int      *Flag3I;
	struct tupla **Tupla2I;
	struct tupla **Tupla3I;
	int      *Flag2IRev;
	int      *Flag3IRev;
	struct tupla **Tupla2IRev;
	struct tupla **Tupla3IRev;
	long n,m;
	double *u;    // u[j]*rcr[j] is the "base" weighted utility of story j
	double *p;    // p[j]*run[j] is the number of story points of story j
	double *pmax; // pmax[i] is the capacity of sprint i, measured in story points
	long *nY;     // nY[j] is the number of stories similar to story j
	long **Y;     // Y[j][k] is the k-th story similar to story j
	long *nUOR;   // nUOR[j] is the number of stories having an OR dependency with story j
	long **UOR;   // UOR[j][k] is the k-th story having an OR dependency with story j
	long *nUAND;  // nUAND[j] is the number of stories having an AND dependency with story j
	long **UAND;  // UAND[j][k] is the k-th story having an AND dependency with story j
	long *FDepA;  // FDepA[j] is 1 if the story j is involved in some dependency (Active)
	long *FDepP;  // FDepP[j] is 1 if the story j is involved in some dependency (Passive)
	long *Cuts;
	int Pred;     // Does it apply precedence inequalities?
	int Cover;    // Does it apply cover inequalities?
	int Lifting;  // Does it apply lifting?
	int KnapSol;  // How solve the KP?
	unsigned long last_cut_checksum;   // Per evitare duplicati
	double last_activity;              // SOLO debug/log
	double last_rhs;                   // SOLO debug/log
	int MaxCutsPerNode;                // Limite per nodo
	int cuts_added_this_node;          // contatore nodo corrente
	Knapsack KP;
};

typedef struct cutinfo CUTINFO, *CUTINFOptr;

// Function implementing a custom cutting plane procedure called by callback
static int CPXPUBLIC mycutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle, int *useraction_p);
static int CPXPUBLIC myNewcutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle, int *useraction_p);
static int CPXPUBLIC myHeucutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle, int *useraction_p);
static int CPXPUBLIC myHeucutcallback_new(CPXCENVptr env, void* cbdata, int wherefrom,
	void* cbhandle, int* useraction_p);
int OR_AND_Inequalities(CUTINFOptr cutinfo);
int OR_AND_InequalitiesHeu(CUTINFOptr cutinfo);
int Cover_Inequalities(CUTINFOptr cutinfo);
int Cover_InequalitiesCplex(CUTINFOptr cutinfo);
int Precedence_Inequalities(CUTINFOptr cutinfo);
int Precedence_Inequalities_2Item(CUTINFOptr cutinfo);
int Precedence_Inequalities_2ItemRev(CUTINFOptr cutinfo);
int Precedence_Inequalities_3Item(CUTINFOptr cutinfo);
int Precedence_Inequalities_3ItemRev(CUTINFOptr cutinfo);
int CalcolaFlag2I(CUTINFOptr cutinfo);
int CalcolaFlag3I(CUTINFOptr cutinfo);
int CalcolaFlag2IRev(CUTINFOptr cutinfo);
int CalcolaFlag3IRev(CUTINFOptr cutinfo);

class CplexObj
{
private:

	int status;

	int CplexObj::SetupKPDP(CUTINFOptr cutinfo);
	int CplexObj::SetupKPMIP(CUTINFOptr cutinfo);

public:

	CPXENVptr env;  // Pointer to Cplex Environment 
	CPXLPptr lp;    // Pointer to LP instance
	CUTINFO cutinfo;  // Data structure used to add cutting planes
	
	// Input
	int minmax;     // Mininum or Maximum (1=minimize; -1=maximize) 
	int ncols;      // Number of columns/variables
	int nrows;      // Number of rows/constraints
	int nz;         // Number of non-zero coefficients into the constraint matrix
	int objsen;     // Minimization/Maximization Flag (+1=Min; -1=Max)
	double *obj;    // Objective Function coefficients
	double *rhs;    // Right Hand Side
	char *sense;    // Constraint Senses
	int *matbeg;    // Pointers to the beginning of each columns of the constraint matrix
	int *matcnt;    // Cardinality of each columns of the constraint matrix
	int *matind;    // Row index associated to each coefficient of the constraint matrix
	double *matval; // Value of each coefficient of the constraint matrix
	double *lb;     // Lower Bound value of each variable
	double *ub;     // Upper Bound value of each variable
	int *indices;
	int *priority;
	int *direction;
	char *xctype;
	double *rngval;

	int nsos;       // Number of SOS (Special Ordered Set) 
	int nsosnz;     // Total number of membres in all SOS added
	char *typesos;  // Array containing SOS type (CPX_TYPE_SOS1='1', CPX_TYPE_SOS2='2')
	int *sosbeg;    // Pointers to the beginning of each SOS
	int *sosind;    // Index associated to the SOS member
	double *soswt;  // Weight associated to the SOS member

	// Output
	double *xt;     // Variable Values
	double objval;  // Objective Function Value

	// Constructor and Destructor
	CplexObj(void);
	~CplexObj(void);

	// MallocCols: Allocate the "column" data structure
	int CplexObj::MallocCols(void);

	// CopyLP
	int CplexObj::CopyLP(void);

	// Set a Lower Bound for the MIP
	int CplexObj::SetLB(double LB);

	// SetMIP
	int CplexObj::SetMIP(double TLim);

	// Sentinel Constraint
	int CplexObj::Sentinel(long n, long m);

	// SolveMIP
	int CplexObj::SolveMIP(double *Zopt, double *Gap, long *Nodes, long *Cuts);

    // Setup Custom Cutting Plane
	int CplexObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax, 
		                       long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP, 
							   int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts);
	int CplexObj::CuttingPlaneHeu(long n, long m, double* u, double* p, double* pmax, 
  	                          long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	                          int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts, long fact);
	int CplexObj::CuttingPlaneHeu_new(long n, long m, double* u, double* p, double* pmax,
		long* nY, long** Y, long* nUOR, long** UOR, long* nUAND, long** UAND, long* FDepA, long* FDepP,
		int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long* Cuts, long fact);
};