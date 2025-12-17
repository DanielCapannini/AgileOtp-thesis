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
#include <windows.h>
#include <processthreadsapi.h>
#include <string>
#include <iomanip>



const bool  WITHOUT_CONSTRAINTS = false; //FLAG TO RUN THE PROBLEM WITHOUT DEPENDENCY CONSTRAINTS

const bool MODIFY_DATA_FOR_TESTING = false; // FLAG TO USE TO USE THE CAPACITY AND SPRINT NUMBER CHANGES

const bool GROUP_SECOND_GOAL_SPRINTS = false; // FLAG TO USE TO GROUP THE SPRITTS OF THE SECOND HALF OF THE PROJECT IN 2 BY 2

int number_of_sprints_to_reduce = 0; // number of sprints to modify for testing
float percentage_of_reduced_story_points_per_modified_sprint = 0.2f; // each modified sprint will have its capacity reduced by 2 story points

long MangiaPath(char *fname);

void main(void)
{
	FILE *fpro;
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

	std::string copy_of_filename;

	if(!MODIFY_DATA_FOR_TESTING && !GROUP_SECOND_GOAL_SPRINTS)
		copy_of_filename = "standard";
	else
		copy_of_filename = "modified";


	if( !WITHOUT_CONSTRAINTS ){
		printf("Running AgileOpt without dependency constraints\n");
		fout1 = fopen("..\\Agile-Heu.out","w");
		fout2 = fopen(("..\\Agile " + copy_of_filename + "-Pro.out").c_str(),"w");
	}
	else {
		printf("Running AgileOpt with dependency constraints\n");
		fout1 = fopen("..\\AgileWithoutDep-Heu.out","w");
		fout2 = fopen(("..\\AgileWithoutDep " + copy_of_filename + "-Pro.out").c_str(),"w");
	}
	printf("\nsono arrivato -1\n");


	// Heuristics
	fprintf(fout1,"\\begin{tabular}{|ll|rrrrr|rrr|rrr|rrrr|rrrr|} \n");
	fprintf(fout1,"\\hline\n");
	fprintf(fout1,"\n\\multicolumn{2}{c|}{   } & ");
	fprintf(fout1,"\n\\multicolumn{5}{c|}{CPLEX} &     ");
	fprintf(fout1,"\n\\multicolumn{3}{c|}{Heuristic} &    ");
	fprintf(fout1,"\n\\multicolumn{3}{c|}{Improved Heuristic} &      ");
	fprintf(fout1,"\n\\multicolumn{4}{c|}{Lagrangian Heuristic} &    ");
	fprintf(fout1,"\n\\multicolumn{4}{c}{Lagrangian Heuristic (coef)} \\\\  \n");
	fprintf(fout1,"\\hline\n");

	fprintf(fout1,"Istanza        &     &             ");
	fprintf(fout1,"  Zopt &     Gap &     Nodes &   NCuts &    Time &   ");
	fprintf(fout1,"  Zheu &      Time &    GHeu &    ");
	fprintf(fout1,"  Zheu &      Time &    GHeu &    ");
	fprintf(fout1,"  Zheu &      Time &    Gap &    GHeu &    ");
	fprintf(fout1,"  Zheu &      Time &    Gap &    GHeu \\\\   \n");
	fflush(fout1);

	// Print instance characteristics
	fprintf(fout2,"Istanza                         ");
	fprintf(fout2," &    &   n &      m &     min(pmax) &     max(pmax) &     avg(pmax) &     nY &    |U^OR| &  |U^AND| \\\\  \n");
	fflush(fout2);
	 
	fscanf(fpro,"%d",&npro);
	printf("\nsono arrivato 0\n");

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
		printf("\nsono arrivato 1\n");

		Prob.ReadData(fname);

		

		j=MangiaPath(fname);
		fprintf(fout1,"%-29s  &     & ",fname+j);
		fprintf(fout2,"%-29s  &     & ",fname+j);
		fprintf(Prob.ferr,"\n%s:\n\n",fname);
		fflush(Prob.ferr);
		printf("\nsono arrivato 2\n");

		// Print instance characteristics
		kount1=0;
		kount2=0;
		kount3=0;
		if(WITHOUT_CONSTRAINTS){
			for (j=0; j<Prob.n; j++)
			{
				Prob.nUOR[j]=0;
				Prob.nUAND[j]=0;
			}
		}
		for (j=0; j<Prob.n; j++)
		{
			if (Prob.nY[j]>1) kount1++;
			if (Prob.nUOR[j]>0) kount2++;
			if (Prob.nUAND[j]>0) kount3++;
		}
		
		if (MODIFY_DATA_FOR_TESTING) {
			if(number_of_sprints_to_reduce > Prob.m) {
				Prob.m = 0;
			}else {
				Prob.m = Prob.m - number_of_sprints_to_reduce;
			}
			for(j = 0; j<Prob.m; j++) {
				Prob.pmax[j] = Prob.pmax[j] * (1.0 - percentage_of_reduced_story_points_per_modified_sprint);
			}
		}

		if (GROUP_SECOND_GOAL_SPRINTS) {
			int half_m = (Prob.m / 2) %2 ==0 ? Prob.m /2 : (Prob.m /2)+1;
			Prob.m = Prob.m - half_m + (half_m /2);
			for(j = 0; j < half_m; j+=2) {
				Prob.pmax[half_m + ((j-1)/2) +1] = Prob.pmax[half_m + j] + Prob.pmax[half_m + j +1];
			}
		}

		double min_pmax = Prob.pmax[0];
		double max_pmax = Prob.pmax[0];
		double sum_pmax = Prob.pmax[0];

		for (int i = 1; i < Prob.m; i++) {
			if (Prob.pmax[i] < min_pmax)
				min_pmax = Prob.pmax[i];
			if (Prob.pmax[i] > max_pmax)
				max_pmax = Prob.pmax[i];
			sum_pmax += Prob.pmax[i];
		}

		double avg_pmax = sum_pmax / Prob.m;

		fprintf(fout2,"%8d & ",Prob.n);
		fprintf(fout2,"%8d & ",Prob.m);
		fprintf(fout2,"%8.1f & ",min_pmax);
		fprintf(fout2,"%8.1f & ",max_pmax);
		fprintf(fout2,"%8.1f & ",avg_pmax);
		fprintf(fout2,"%8d & ",kount1);
		fprintf(fout2,"%8d & ",kount2);
		fprintf(fout2,"%8d \\\\ \n",kount3);
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
		dt = (double)(t2 - t1) / CLOCKS_PER_SEC;

		zlp = Prob.Zopt;


		fprintf(fout1,"%10.1lf & ",Prob.Zopt);
		fprintf(fout1,"%7.2lf & ",Prob.Gap);
		fprintf(fout1,"%8d & ",Prob.Nodes);
		fprintf(fout1,"%8d & ",Prob.Cuts);
		fprintf(fout1,"%8.2lf & ",dt);
		fflush(fout1);

		// Heuristic procedure
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;


		t1 = clock();
		Prob.OptimizeHeu();
		t2 = clock();
		dt = (double)(t2 - t1) / CLOCKS_PER_SEC;

		fprintf(fout1,"%10.1lf & ",Prob.Zheu);
		fprintf(fout1,"%8.2lf & ",dt);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf & ",gap);
		fflush(fout1);
		zheu = Prob.Zheu;


		Prob.cfg.TimeLimit = TLimHEU;
		Prob.cfg.Sentinel = 0;

		t1 = clock();
		Prob.OptimizeHeu_Improved();
		t2 = clock();
		dt = (double)(t2 - t1) / CLOCKS_PER_SEC;

		fprintf(fout1, "%10.1lf & ", Prob.Zheu);
		fprintf(fout1, "%8.2lf & ", dt);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf & ",gap);
		fflush(fout1);
		zheu = Prob.Zheu;

		// Lagrangian Heuristic procedure
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;
		Prob.cfg.HeuBest=1.0;


		t1 = clock();
		Prob.OptimizeLagrHeu();
		t2 = clock();
		dt = (double)(t2 - t1) / CLOCKS_PER_SEC;


		fprintf(fout1,"%10.1lf & ",Prob.Zheu);
		fprintf(fout1,"%8.2lf & ",dt);
		gap = 100.*((Prob.Zlagr-Prob.Zheu)/Prob.Zheu);
		fprintf(fout1,"%7.2lf & ",gap);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf & ",gap);
		fflush(fout1);
		zheu = Prob.Zheu;

		// Lagrangian Heuristic procedure with improvement even with HeuBest*Zub solution
		Prob.cfg.TimeLimit=TLimHEU;
		Prob.cfg.Sentinel=0;
		Prob.cfg.HeuBest=0.995; // Best=0.995


		t1 = clock();
		Prob.OptimizeLagrHeu();
		t2 = clock();
		dt = (double)(t2 - t1) / CLOCKS_PER_SEC;

		fprintf(fout1,"%10.1lf & ",Prob.Zheu);
		fprintf(fout1,"%8.2lf & ",dt);
		gap = 100.*((Prob.Zlagr-Prob.Zheu)/Prob.Zheu);
		fprintf(fout1,"%7.2lf & ",gap);

		gap = 100.*((zlp-Prob.Zheu)/zlp+0.00001);
		fprintf(fout1,"%7.2lf \\\\ \n",gap);
		fflush(fout1);
		zheu = Prob.Zheu;
	}

	fclose(fpro);
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