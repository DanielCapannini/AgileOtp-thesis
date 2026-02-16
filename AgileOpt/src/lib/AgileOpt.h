//-----------------------------------------------------------------------------
//  File: AgileOpt.h                                                         
//                                                                           
//  Project: Optimize the Agile Project Management Schedule (header)                 
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 22.08.2011                                                  
//-----------------------------------------------------------------------------
//

#define  Prec    0.00001
#define  LPrec   0.99999
#define  Inf     10000000000.
#define  iprinx  0

struct config
{
	double TimeLimit;  // Total time limit
	int Sentinel;      // Sentinel=1 if the Sentinel Constraints must be applied
	int Pred;          // Pred=1 if the Dominance Inequalities must be applied
	int Cover;         // Cover=1 if the Cover Inequalities must be applied
	int Lifting;       // Lifting=1 if the Lifted Cover Inequalities must be applied
	int KnapSol;       // Knapsack=0 use Dynamic Programming; Knapsack=1 use Cplex    
	int AddLB;         // AddLB=0 do not add a lower bound to the MIP; AddLB=1 add!
	int MaxCuts;       // Maximum number of user cuts
	double HeuBest;       // Percentage of Best Zub to consider to apply the strong greedy heuristic
};

class AgileOpt
{
private:
	long GetList(FILE *fin, long *nlist, long *list);
	int AgileOpt::Reduction(void);
	int AgileOpt::Feasibility(double *xheu, double *yheu, double *zheu, int *bol);
	int AgileOpt::ComputePenPrice(double *prpen, double *Lambda, double *LambdaOR, double *LambdaAND, double *Konst);
	int AgileOpt::ComputePenPriceDP(double *prpen,double *apen, double *Lambda, double *LambdaOR, double *LambdaAND, 
		                            double *LambdaY1, double *LambdaY2, double *Konst);
	int AgileOpt::SolveLR(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double Konst, double *zlr, double *xlr, double *ylr);
	int AgileOpt::SolveLRDP(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double Konst, double *zlr, double *xlr, double *ylr);
	double AgileOpt::TryMove(double *xheu, double *yheu, long j1, long i1, long i);
	double AgileOpt::Move(double *xheu, double *yheu, long j1, long i1, long i);
	int AgileOpt::HeuSwaps(double *zheu, double *xheu, double *yheu);
	int AgileOpt::HeuSwapsOld(double *zheu, double *xheu, double *yheu);
	int AgileOpt::LagrHeu(double *urpen, double *apen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheu, double *yheu, int select_method);
	int AgileOpt::LagrHeuDP(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheu, double *yheu);
	int AgileOpt::LagrHeuDPBack(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheub, double *yheub, int select_method);
	int AgileOpt::LagrHeuDP_Old(double *urpen, long nYMap, long *YMap, long *Flag, double *zheu, double *xheu);
	int AgileOpt::UpdatePen(double alpha, double zheu, double zlr, double *xlr, double *Lambda, double *LambdaOR, double *LambdaAND, 
		                    double *subgr, double *subgrOR, double *subgrAND);
	int AgileOpt::UpdatePenDP(double alpha, double zheu, double zlr, double *xlr, double *ylr, double *Lambda, double *LambdaOR, double *LambdaAND, 
	                      double *LambdaY1, double *LambdaY2, double *subgr, double *subgrOR, double *subgrAND, double *subgrY1, double *subgrY2);

public:

	FILE *ferr;   // Log Error!

	long n;       // Number of User Stories
	long m;       // Number of Sprints
	long original_m; // Original number of Sprints before any modification
	double *a;    // a[j] is the affinity of story j
	double *u;    // u[j] is the utility of story j
	double *ur;   // u[j]*rcr[j] is the "base" weighted utility of story j
	double *p;    // p[j] is the number of story points of story j
	double *pr;   // p[j]*run[j] is the weight of story j
	double *pmax; // pmax[i] is the capacity of sprint i, measured in story points
	double *rcr;  // rcr[j] is the criticality risk of story j
	double *run;  // run[j] is the uncertainty risk of story j
	long *nY;     // nY[j] is the number of stories similar to story j
	long **Y;     // Y[j][k] is the k-th story similar to story j
	long *nD;     // nD[j] is the number of stories the story j depends on
	long **D;     // D[j][k] is the k-th story the story j depends on
	long *nUOR;   // nUOR[j] is the number of stories having an OR dependency with story j
	long **UOR;   // UOR[j][k] is the k-th story having an OR dependency with story j
	long *nUAND;  // nUAND[j] is the number of stories having an AND dependency with story j
	long **UAND;  // UAND[j][k] is the k-th story having an AND dependency with story j
	long *FDepA;  // FDepA[j] is 1 if the story j is involved in some dependency (active)
	long *FDepP;  // FDepP[j] is 1 if the story j is involved in some dependency (passive)

	struct config cfg;  // Configuration of the exact algorithm
	double Zopt;
	double Zheu;
	double Zlagr;
	double Gap;
	long Nodes;
	long Cuts;
	double NSprintsUsed;
	double PercentageUtilization;
	double DeviationRisk;
	double DeviationUncertaintyRisk;
	double UncertaintyMin;
	double UncertaintyMax;
	double RiskMin;
	double RiskMax;
	double HalfUtilitiSprint;

	// Constructor and Destructorr
	AgileOpt(void);
	~AgileOpt(void);

	// Read the instance from Text File "name"
	int ReadData(const char *name); 
	int ReadDataWithoutConstraints(const char *name);

	// Optimize the Agile Schedule
	int AgileOpt::Optimize(int select_method, char* fname);
	int AgileOpt::Optimize2(int select_method, char* fname);
	int AgileOpt::OptimizeHeu(int select_method, char* fname);
	int AgileOpt::OptimizeHeu_Improved(int select_method, char* fname);
	int AgileOpt::OptimizeLagrHeu(int select_method, char* fname);
};


