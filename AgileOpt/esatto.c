//-----------------------------------------------------------------------------
//  File: esatto.c
//
//  Program to assign customers to areas 
//
//  Authors: M.A.Boschetti, V.Maniezzo          
//
//  Last update: 24.11.2010 
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <tchar.h>
#include <time.h>
#include "heumed.h"
#include "cplex.h"

#define INF  1000000000

// ----------------------------------------------------------------------
int esatto(int BMIP, struct problema *pro, long M, long *MED, long *PMED, long **PRO, long **WEIT, long *CAPAC,
		   long *ZZ, long *XSTAR, long *XSTARV)
{
	CPXENVptr env=NULL;
	CPXLPptr lp=NULL;

	int status;
    char *probname = "lp";
	int ncols;
	int nrows;
	int objsen;
	double objval;
	double *obj;
	double *rhs;
	char *sense;
	int *matbeg;
	int *matcnt;
	int *matind;
	double *matval;
	double *lb;
	double *ub;
	double *xt,XX,XXt;
	int *indices;
	char *xctype;
	double *rngval=NULL;

	// Auxiliary variables
	long kc,kcc,knz,kcp;
	long nvar1,nvar2,nvar3,nvar4;
	long ncon1,ncon2,ncon3,ncon4,ncon5,ncon6;
	long ncon1a,ncon1b,ncon1c;
	int param;
	int ivalue;
	double rvalue;
	double distx,disty;
	double tdep;
	long i,j,k,i1,j1,idep;
	long ii;
	long nk;
	long qtot,min,max,avg;
	//int BMIP=1;
	//int BMIP=0;
    double tt;
    clock_t t1,t2;

	// Cutting plane
//	CUTINFO cutinfo;

	t1 = clock();

	// Open Cplex Enviroment
	if (env==NULL)
	{
		env = CPXopenCPLEX(&status);
		if (status)
		{
			printf("\nError... CPXopenCPLEX...\n");
			return 1;
		}

		lp = CPXcreateprob(env,&status,probname);
		if (status)
		{
			printf("\nError... CPXcreateprob...\n");
			return 2;
		}
	}

	param=CPX_PARAM_DATACHECK;
	ivalue=1;
	//ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_SCRIND;
	ivalue=1;
	// ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_MIPDISPLAY;
	ivalue=2;
	// ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_SIMDISPLAY;
	//ivalue=2;
	ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	// Setup data

	// Variables
	nk = pro->vehi.ntipi;
	nvar1 = pro->n*M*nk;   // variables x_ijk 
	nvar2 = M*nk;          // variables y_jk 
	ncols = nvar1+nvar2;
	//ncols = nvar1+nvar2+nvar3;

	// Constraints
	ncon1a = M*nk;      // Capacity constraints: peso 
	ncon1b = M*nk;      // Capacity constraints: volume 
	ncon1c = M*nk;      // Capacity constraints: tempo
	ncon1=ncon1a+ncon1b+ncon1c;
	ncon2 = pro->n;    // Assignement of each customer to exactly one area
	ncon3 = nk;        // Maximum number of vehicles to use 
	ncon4 = M;         // Assignment of each area to exactly one vehicle
	ncon5 = 2;         // Minimum and Maximum number of areas
	nrows = ncon1+ncon2+ncon3+ncon4+ncon5;

	// Minimization
	objsen = +1;

	obj = (double*)malloc(ncols*sizeof(double));
	rhs = (double*)malloc(nrows*sizeof(double));
	sense = (char*)malloc(nrows*sizeof(char));
	matbeg = (int*)malloc((ncols+4)*sizeof(int));
	matcnt = (int*)malloc((ncols+4)*sizeof(int));
	matind = (int*)malloc(6*ncols*sizeof(int));          // DA VERIFICARE!
	matval = (double*)malloc(6*ncols*sizeof(double));       // DA VERIFICARE!
	lb = (double*)malloc(ncols*sizeof(double));
	ub = (double*)malloc(ncols*sizeof(double));
	xt = (double*)malloc(ncols*sizeof(double));
	indices = (int*)malloc(ncols*sizeof(int));
	xctype = (char*)malloc(ncols*sizeof(char));

	// Variables x_ijk
	kc=0;
	knz=0;
	for (ii=1;ii<=pro->n;ii++)
	{
		i=pro->cust.sind[ii];
		for (j=1;j<=M;j++)
		{
			for (k=1;k<=nk;k++)
			{
				// Objective function
				//obj[kc]=(double)PRO[j][i];

				obj[kc]=(double)Dist(pro,i,MED[j]);
				//if (obj[kc]<-0.001)
				//	kc=kc;

				// Column
				matbeg[kc]=knz;
				matcnt[kc]=0;

				// Capacity constraints: peso
				matcnt[kc]++;
				matind[knz]=(j-1)*nk+k-1;
				//matval[knz]=(double)WEIT[j][i];
				if (pro->cust.peso[i]<=pro->vehi.peso[k])
					matval[knz]=(double)pro->cust.peso[i];
				else
					matval[knz]=(double)pro->vehi.peso[k]; // Errore nei dati che segnaleremo dopo
				knz++;

				// Capacity constraints: volume
				matcnt[kc]++;
				matind[knz]=ncon1a+(j-1)*nk+k-1;
				//matval[knz]=(double)WEIT[j][i];
				if (pro->cust.volume[i]<=pro->vehi.volume[k])
					matval[knz]=(double)pro->cust.volume[i];
				else
					matval[knz]=(double)pro->vehi.volume[k]; // Errore nei dati che segnaleremo dopo
				knz++;

				// Capacity constraints: tempo
				matcnt[kc]++;
				matind[knz]=ncon1a+ncon1b+(j-1)*nk+k-1;
				//matval[knz]=(double)WEIT[j][i];
				if (pro->cust.tempo[i]<=pro->vehi.tempo[k])
					matval[knz]=(double)pro->cust.tempo[i];
				else
					matval[knz]=(double)pro->vehi.tempo[k]; // Errore nei dati che segnaleremo dopo
				knz++;

				// Assignment constraints 
				matcnt[kc]++;
				matind[knz]=ncon1+ii-1;
				matval[knz]=+1.0;
				knz++;

				// Lower and upper bound
				lb[kc]=0.;
				if ((pro->cust.mezzo[i]>0)&&(pro->cust.mezzo[i]!=k))
					ub[kc]=0.;
				else
					ub[kc]=1.;
				indices[kc]=kc;
				//xctype[kc]='B';
				xctype[kc]='C';

				kc++;
			}
		}
	}

	// Variables y_jk
	for (j=1;j<=M;j++)
	{
		for (k=1;k<=nk;k++)
		{
			// Objective function
			if ((pro->flagdep)&&(pro->flagobj==1))
				obj[kc] = (double)DistDep(pro,MED[j]);
			 else
				obj[kc]=0.0;   

			// Column
			matbeg[kc]=knz;
			matcnt[kc]=0;

			// Capacity constranits: peso
			matcnt[kc]++;
			matind[knz]=(j-1)*nk+k-1;
			//matval[knz]=-(double)CAPAC[j];
			matval[knz]=-(double)pro->vehi.peso[k];
			knz++;

			// Capacity constranits: volume
			matcnt[kc]++;
			matind[knz]=ncon1a+(j-1)*nk+k-1;
			//matval[knz]=-(double)CAPAC[j];
			matval[knz]=-(double)pro->vehi.volume[k];
			knz++;

			// Capacity constranits: tempo
			matcnt[kc]++;
			matind[knz]=ncon1a+ncon1b+(j-1)*nk+k-1;
			//matval[knz]=-(double)CAPAC[j];
			if (pro->flagtdep)
			{
				j1=pro->cust.ptr[MED[j]];
				idep=pro->dep.ptr;
				tdep = (double)min(pro->mat.time[idep][j1],pro->mat.time[j1][idep]);
				matval[knz]=-(double)pro->vehi.tempo[k]+tdep;
			}
			else
				matval[knz]=-(double)pro->vehi.tempo[k];
			knz++;

			// Fleet constraints: maximum
			matcnt[kc]++;
			matind[knz]=ncon1+ncon2+k-1;
			matval[knz]=1.;
			knz++;

			// Assignment constraints
			matcnt[kc]++;
			matind[knz]=ncon1+ncon2+ncon3+j-1;
			matval[knz]=1.;
			knz++;

			// Minimum number of areas
			matcnt[kc]++;
			matind[knz]=ncon1+ncon2+ncon3+ncon4;
			matval[knz]=1.;
			knz++;

			// Maximum number of areas
			matcnt[kc]++;
			matind[knz]=ncon1+ncon2+ncon3+ncon4+1;
			matval[knz]=1.;
			knz++;

			// Lower and upper bound
			lb[kc]=0.0;
			ub[kc]=1.0;

			indices[kc]=kc;
			xctype[kc]='B';

			kc++;
		}
	}

	matbeg[kc]=knz;
	ncols=kc;

	// Right Hand Side: Capacity Constraints (3x)
	ii=0;
	for (j=1;j<=ncon1;j++)
	{
		rhs[ii]=0.0;
		sense[ii]='L';
		ii++;
	}

	// Right Hand Side: Assignment Constraints
	for (i=1;i<=pro->n;i++)
	{
		rhs[ii]=+1.0;
		sense[ii]='E';
		ii++;
	}

	// Right Hand Side: Fleet Constraints
	for (k=1;k<=nk;k++)
	{
		rhs[ii]=(double)pro->vehi.num[k];
		sense[ii]='L';
		ii++;
	}

	// Right Hand Side: Assignment Constraints
	for (j=1;j<=M;j++)
	{
		rhs[ii]=1.0;
		sense[ii]='L';
		ii++;
	}

	// Right Hand Side: Minimum number of areas
	rhs[ii]=pro->nummed1;
	sense[ii]='G';
	ii++;

	// Right Hand Side: Maximum number of areas
	rhs[ii]=pro->nummed2;
	sense[ii]='L';
	ii++;

	rngval=NULL;

	// Load the data into Cplex
    status = CPXcopylp(env,lp,ncols,nrows,objsen,obj,rhs,sense,matbeg,matcnt,matind,matval,lb,ub,rngval);

	param=CPX_PARAM_TILIM;
	if (BMIP)
		rvalue=1000.;
	else
		rvalue=100.;
	status = CPXsetdblparam (env, param, rvalue);
	if (status>0) return status;

	param=CPX_PARAM_THREADS;
	ivalue=1;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

// Callback for cut generation

  // /* Set parameters */

  // /* Assure linear mappings between the presolved and original
  //    models */

  // status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
  // if (status) goto TERMINATE;


  // /* Turn on traditional search for use with control callbacks */

  // status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
  // if (status)  goto TERMINATE;

  // /* Let MIP callbacks work on the original model */

  // status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
  // if ( status ) goto TERMINATE;

  // cutinfo.pro=pro;
  // cutinfo.Flag=1;
  // if (BolY==0)
  // {
		//cutinfo.start=nvar1+nvar2;
		//cutinfo.end=nvar1+nvar2+nvar3;
  // }
  // else
  // {
		//cutinfo.start=0;
		//cutinfo.end=nvar1;
  // }
  // cutinfo.lp = lp;
  // cutinfo.lps = NULL;
  // cutinfo.MatX = MatX;
  // cutinfo.MatY = MatY;
  // cutinfo.MatXi = MatXi;
  // cutinfo.numcols = CPXgetnumcols (env, lp);
  // cutinfo.x = (double *) malloc (cutinfo.numcols * sizeof (double));
  // cutinfo.beg = (int *) malloc (2 * sizeof (int));
  // cutinfo.ind = (int *) malloc (cutinfo.numcols * sizeof (int)); 
  // cutinfo.val = (double *) malloc (cutinfo.numcols * sizeof (double));
  // cutinfo.rhs = (double *) malloc (1 * sizeof (double));
  // cutinfo.indices = (int *) malloc (cutinfo.numcols * sizeof (int));
  // cutinfo.values = (double *) malloc (cutinfo.numcols * sizeof (double));
  // if ((cutinfo.x == NULL)||(cutinfo.beg == NULL)||(cutinfo.ind == NULL)||
	 //  (cutinfo.val == NULL)||(cutinfo.rhs == NULL)||(cutinfo.indices == NULL)||(cutinfo.values == NULL)) 
  // {
  //    fprintf (stderr, "No memory for solution values.\n");
  //    goto TERMINATE;
  // }

  // /* Create user cuts for noswot problem */

  // if (BolY==0)
	 //  status = BuildSepModel (env, &cutinfo);
  // else
	 //  status = BuildSepModelY(env, &cutinfo);
  // if ( status )  goto TERMINATE;

  // /* Set up to use MIP callback */

  // if (BolY==0)
		//status = CPXsetcutcallbackfunc (env, mycutcallback, &cutinfo);
  // else
		//status = CPXsetcutcallbackfunc(env, mycutcallbackY, &cutinfo);

  // if ( status )  goto TERMINATE;

// End Callback for cut generation

	if (BMIP)
	{
		// Solve the MIP
		printf("\n--------\nSolve the MIP\n\n");
		status = CPXchgprobtype(env,lp,CPXPROB_MILP);
		status = CPXchgctype(env,lp,ncols,indices,xctype);
		//status = CPXwriteprob(env,lp,"myprob.lp","LP");

		// Solve the problem by Cplex
		status = CPXmipopt(env,lp);
		//printf("MIPOpt Status: %d\n\n",status);

		status = CPXgetstat(env,lp);
		if ((status!=CPXMIP_OPTIMAL)&&(status!=CPXMIP_OPTIMAL_TOL)&&(status!=CPXMIP_FEASIBLE)&&
			(status!=CPXMIP_TIME_LIM_FEAS))
		{
			printf("Non Ammissibile...\n",status);
			*ZZ=Infi;
			goto TERMINATE;
		}

		//printf("MIP Status: %d\n\n",status);
		//fflush(LogFile);

		// Read the solution
//		status = CPXgetmipobjval(env,lp,&objval);
		status = CPXgetobjval(env,lp,&objval);
		*ZZ=(long)(objval);
//		status = CPXgetmipx(env,lp,xt,0,ncols-1);
 		status = CPXgetx(env,lp,xt,0,ncols-1);
		// status = CPXgetpi(env,lp,u,begin,end);
// 		Display_Solution3(pro,ncols,xt,MatY,MatX,MatXi);

		kc=0;
		knz=0;
		*ZZ=0;
		if ((pro->flagdep)&&(pro->flagobj==1))
		{
			for (j=1;j<=M;j++)
			{
				(*ZZ) += DistDep(pro,MED[j]);
			}
		 }

		for (i=1;i<=pro->n;i++)
		{
			XSTAR[i]=-1;
			XSTARV[i]=-1;
			if (pro->flaground==0)
				XX=xt[kc];
			else
				XX=+Infr;
			for (j=1;j<=M;j++)
			{
				for (k=1;k<=nk;k++)
				{
					if ((xt[kc]>0.001)&&(xt[kc]<0.999))
						knz++;

					if (pro->flaground==0)
					{
						if (xt[kc]>XX)
						{
							XSTAR[i]=j;
							XSTARV[i]=k;
							XX=xt[kc];
						}
					}
					else
					{
						ii=pro->cust.sind[i];
						XXt = Dist(pro,ii,MED[j]);
						XXt = (1.-xt[kc])*XXt;
						if ((xt[kc]>0.001)&&(XXt<XX))
						{
							XSTAR[i]=j;
							XSTARV[i]=k;
							XX=XXt;
						}
					}
					kc++;
				}
			}

			if (XSTAR[i]<0)
			{
				printf("Soluzione errata...\n",status);
				*ZZ=Infi;
				goto TERMINATE;
			}

			ii=pro->cust.sind[i];
			(*ZZ) += Dist(pro,ii,MED[XSTAR[i]]);

		}

		t2 = clock();
		tt = ((double)(t2-t1))/CLK_TCK;///CLOCKS_PER_SEC;

		printf("MIP obj: %12.2lf / %d\n",objval,(*ZZ));
		printf("--> Time: %5.2lf \t n.var.fraz.: %d\n",tt,knz);

	}
	else
	{
		// Solve the LP
		printf("\n--------\nSolve the LP\n\n");
		// status = CPXwriteprob(env,lp,"myprob.lp","LP");
		// param=CPX_PARAM_SCRIND;
		// ivalue=1;
		// status = CPXsetintparam (env, param, ivalue);
		// if (status>0) return status;

		// Solve the problem by Cplex
		status = CPXdualopt(env,lp);
		printf("DualOpt Status: %d \t",status);

		status = CPXgetstat(env,lp);
		printf("LP Status: %d\t",status);
		//fflush(LogFile);
		if (status!=CPX_STAT_OPTIMAL)
		{
			printf("Non Ammissibile...\n",status);
			*ZZ=Infi;
			goto TERMINATE;
		}

		// Read the solution
		status = CPXgetobjval(env,lp,&objval);
		*ZZ=(long)(objval);
		status = CPXgetx(env,lp,xt,0,ncols-1);
// 		Display_Solution3(pro,ncols,xt,MatY,MatX,MatXi);

		kc=0;
		knz=0;
		*ZZ=0;
		if ((pro->flagdep)&&(pro->flagobj==1))
		{
			for (j=1;j<=M;j++)
			{
				(*ZZ) += DistDep(pro,MED[j]);
			}
		 }

		for (i=1;i<=pro->n;i++)
		{
			XSTAR[i]=1;
			XSTARV[i]=1;
			if (pro->flaground==0)
				XX=xt[kc];
			else
				XX=+Infr;
			for (j=1;j<=M;j++)
			{
				for (k=1;k<=nk;k++)
				{
					if ((xt[kc]>0.001)&&(xt[kc]<0.999))
						knz++;

					if (pro->flaground==0)
					{
						if (xt[kc]>XX)
						{
							XSTAR[i]=j;
							XSTARV[i]=k;
							XX=xt[kc];
						}
					}
					else
					{
						ii=pro->cust.sind[i];
						XXt = Dist(pro,ii,MED[j]);
						XXt=(1.-xt[kc])*XXt;
						if ((xt[kc]>0.001)&&(XXt<XX))
						{
							XSTAR[i]=j;
							XSTARV[i]=k;
							XX=XXt;
						}
					}
					kc++;
				}
			}

			ii=pro->cust.sind[i];
			(*ZZ) += Dist(pro,ii,MED[XSTAR[i]]);

		}

		t2 = clock();
		tt = ((double)(t2-t1))/CLK_TCK;///CLOCKS_PER_SEC;

		printf("LP obj: %12.2lf / %d\n",objval,(*ZZ));
		printf("--> Time: %5.2lf \t n.var.fraz.: %d\n",tt,knz);

	}

/*	// Sistema XSTAR
	for (j=1;j<=M;j++)
		PMED[j]=0;

	for (i=1;i<=pro->n;i++)
		PMED[XSTAR[i]]=1;

	ii=1;
	for (j=1;j<=M;j++)
		if (PMED[j]>0)
		{
			PMED[j]=ii;
			ii++;
		}

	for (i=1;i<=pro->n;i++)
		XSTAR[i]=PMED[XSTAR[i]];
*/
	//z=objval;

TERMINATE:

	// Close the Cplex Enviroment
	status = CPXfreeprob(env, &lp);
	if (status)
	{
		printf("\nError... CPXfreeprob...\n");
		return 1;
	}

	status = CPXcloseCPLEX(&env);
	if (status)
	{
		printf("\nError... CPXcloseCPLEX...\n");
		return 2;
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
	free(xctype);

	//delete [] cutinfo.x;
	//delete [] cutinfo.beg;
	//delete [] cutinfo.ind; 
	//delete [] cutinfo.val;
	//delete [] cutinfo.rhs;
	//delete [] cutinfo.indices;
	//delete [] cutinfo.values;

	return 0;
}