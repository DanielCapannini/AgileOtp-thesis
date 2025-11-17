//-----------------------------------------------------------------------------
//  File: main.cpp                                                           
//                                                                           
//  Project: Optimize the Agile Project Management Schedule                  
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 28.08.2012                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <time.h>
#include "lib\AgileOpt.h"

long MangiaPath(char *fname);

void main(void)
{
	FILE *fpro;
	FILE *fout;
	FILE *fout1;
	FILE *fout2;
	AgileOpt Prob;
	
	long npro;
	long i,j;
	long kount1,kount2,kount3;
	long t1,t2;
	double dt;
	double zheu,gap,zlp;
	char fname[255];
	double TLimMIP;
	double TLimMIPb;
	double TLimHEU;

	fpro = fopen("..\\Agile.pro","r");
	fout = fopen("..\\Agile-LP.out","w");
	fout1 = fopen("..\\Agile-Heu.out","w");
	fout2 = fopen("..\\Agile-Pro.out","w");

	// LP
	fprintf(fout,"\n                              ");
	fprintf(fout,"                   Cplex                 ");
	fprintf(fout,"             Cplex + LagrHeu             ");
	fprintf(fout,"                 Sentinel                ");
	fprintf(fout,"                   Pred                  ");
	fprintf(fout,"                   Cover                 ");
	fprintf(fout,"               Lifted Cover            \n");

	fprintf(fout,"Istanza                         ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout,"  Zopt    Gap    Nodes   NCuts    Time   \n");
	fflush(fout);

	// Heuristics
	fprintf(fout1,"\n                              ");
	fprintf(fout1,"                   Cplex                 ");
	fprintf(fout1,"          Heuristic      ");
	fprintf(fout1,"       Lagrangean Heuristic      ");
	fprintf(fout1,"   Lagrangean Heuristic (coef)   \n");

	fprintf(fout1,"Istanza                         ");
	fprintf(fout1,"  Zopt    Gap    Nodes   NCuts    Time   ");
	fprintf(fout1,"  Zheu     Time   GHeu   ");
	fprintf(fout1,"  Zheu     Time   Gap    GHeu   ");
	fprintf(fout1,"  Zheu     Time   Gap    GHeu   \n");
	fflush(fout1);

	// Print instance characteristics
	fprintf(fout2,"Istanza                         ");
	fprintf(fout2,"    n       m      nY    |U^OR|  |U^AND|   \n");
	fflush(fout2);
	 
	fscanf(fpro,"%d",&npro);

	// Set Time Limit
	Prob.cfg.TimeLimit=10.;
	Prob.cfg.KnapSol=0;
	Prob.cfg.MaxCuts=50000;
	TLimMIP = 600.;
	TLimMIPb = 600.;
	TLimHEU = 600.;

	for (i=0; i<npro; i++)
	{
		fscanf(fpro,"%s",fname);

		Prob.ReadData(fname);
		j=MangiaPath(fname);
		fprintf(fout,"%-29s",fname+j);
		fprintf(fout1,"%-29s",fname+j);
		fprintf(fout2,"%-29s",fname+j);
		fprintf(Prob.ferr,"\n%s:\n\n",fname);
		fflush(Prob.ferr);

		// Print instance characteristics
		kount1=0;
		kount2=0;
		kount3=0;
		for (j=0; j<Prob.n; j++)
		{
			if (Prob.nY[j]>1) kount1++;
			if (Prob.nUOR[j]>0) kount2++;
			if (Prob.nUAND[j]>0) kount3++;
		}
		fprintf(fout2,"%8d",Prob.n);
		fprintf(fout2,"%8d",Prob.m);
		fprintf(fout2,"%8d",kount1);
		fprintf(fout2,"%8d",kount2);
		fprintf(fout2,"%8d\n",kount3);
		fflush(fout2);

		// Cplex
		Prob.cfg.TimeLimit=TLimMIP;
		Prob.cfg.Sentinel=0;
		Prob.cfg.Pred=0;
		Prob.cfg.Cover=0;
		Prob.cfg.Lifting=0;
		Prob.cfg.AddLB=0;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;
		zlp = Prob.Zopt;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		fprintf(fout1,"%10.1lf",Prob.Zopt);
		fprintf(fout1,"%7.2lf",Prob.Gap);
		fprintf(fout1,"%8d",Prob.Nodes);
		fprintf(fout1,"%8d",Prob.Cuts);
		fprintf(fout1,"%8.2lf",dt);
		fflush(fout1);

		// Heuristic procedure
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;

		t1 = clock();
		Prob.OptimizeHeu();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout1,"%10.1lf",Prob.Zheu);
		fprintf(fout1,"%8.2lf",dt);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf",gap);
		fflush(fout1);
		zheu = Prob.Zheu;

		// Lagrangian Heuristic procedure
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;
		Prob.cfg.HeuBest=1.0;

		t1 = clock();
		Prob.OptimizeLagrHeu();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout1,"%10.1lf",Prob.Zheu);
		fprintf(fout1,"%8.2lf",dt);
		gap = 100.*((Prob.Zlagr-Prob.Zheu)/Prob.Zheu);
		fprintf(fout1,"%7.2lf",gap);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf",gap);
		fflush(fout1);
		zheu = Prob.Zheu;

		// Lagrangian Heuristic procedure with improvement even with HeuBest*Zub solution
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;
		Prob.cfg.HeuBest=0.995; // Best=0.995

		t1 = clock();
		Prob.OptimizeLagrHeu();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout1,"%10.1lf",Prob.Zheu);
		fprintf(fout1,"%8.2lf",dt);
		gap = 100.*((Prob.Zlagr-Prob.Zheu)/Prob.Zheu);
		fprintf(fout1,"%7.2lf",gap);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf\n",gap);
		fflush(fout1);
		zheu = Prob.Zheu;

		// Cplex + LagrHeu
		Prob.cfg.TimeLimit=TLimMIPb;
		Prob.cfg.Sentinel=0;
		Prob.cfg.Pred=0;
		Prob.cfg.Cover=0;
		Prob.cfg.Lifting=0;
		Prob.cfg.AddLB=1;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		// Sentinel
		Prob.cfg.TimeLimit=TLimMIPb;
		Prob.cfg.Sentinel=1;
		Prob.cfg.Pred=0;
		Prob.cfg.Cover=0;
		Prob.cfg.Lifting=0;
		Prob.cfg.AddLB=0;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		// Pred
		Prob.cfg.TimeLimit=TLimMIPb;
		Prob.cfg.Sentinel=0;
		Prob.cfg.Pred=1;
		Prob.cfg.Cover=0;
		Prob.cfg.Lifting=0;
		Prob.cfg.AddLB=0;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		// Cover
		Prob.cfg.TimeLimit=TLimMIPb;
		Prob.cfg.Sentinel=0;
		Prob.cfg.Pred=0;
		Prob.cfg.Cover=1;
		Prob.cfg.Lifting=1;
		Prob.cfg.AddLB=0;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		// Lifted Cover
		Prob.cfg.TimeLimit=TLimMIPb;
		Prob.cfg.Sentinel=0;
		Prob.cfg.Pred=1;
		Prob.cfg.Cover=1;
		Prob.cfg.Lifting=1;
		Prob.cfg.AddLB=0;

		t1 = clock();
		Prob.Optimize();
		t2 = clock();
		dt = (double)(t2-t1)/CLOCKS_PER_SEC;

		fprintf(fout,"%10.1lf",Prob.Zopt);
		fprintf(fout,"%7.2lf",Prob.Gap);
		fprintf(fout,"%8d",Prob.Nodes);
		fprintf(fout,"%8d",Prob.Cuts);
		fprintf(fout,"%8.2lf",dt);
		fflush(fout);

		fprintf(fout,"\n");

	}

	fclose(fpro);
	fclose(fout);
	fclose(fout1);
	fclose(fout2);
}

long MangiaPath(char *fname)
{
	long i,j;

	j=0;
	for (i=0; fname[i]!='\0'; i++)
		if (fname[i]==92)
			j=i+1;

	return j;
}