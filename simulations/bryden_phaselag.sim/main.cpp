/******* Phase Lagged Oscillation ********
	Author: Jordan Perr-Sauer
	Summer 2011

	Evolves phase lagged oscillators
as proposed in Bryden's paper.

*****************************************/


#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <ctime>
#include "TSearch.h"
#include "random.h"
#include "VectorMatrix.h"
using namespace std;


/**** Simulation constants ****/
#define MAXTIME 300
#define PRINTSTEP 1
#define SEGMENTS 5

/**** Evolution constants ****/
#define GENERATIONS 150
#define POPULATIONSIZE 60


/**** Convenience macro, return seed as string ****/
#define SEEDSTR ltos(XG_SEED)
string inline ltos(long x){
	char a[50];
	sprintf(a, "%ld", x);
	return a;
}


/**** Convenience function for sigmoid ****/
double inline sigmoid(double x) {
	return 1/(1+exp(-x));
}
#define LIMIT(x) (2*sigmoid(x)-1)


/**** Globals and defaults ****/
// Simulation timestep
double TIMESTEP = 0.01;
// First argument to program. Used as seed and in file names.
long XG_SEED;
// Some internal globals
bool output;
ofstream ofile;

// Accumulators for each segment.
typedef struct {
	double Vvb;
	double Vdb;
	double theta[3];
} Segment;


double simulate(TVector<double> &t, RandomState &rs)
{

	// Create array segment accumulators.
	Segment seg[SEGMENTS];

	// Structs dont get initialized, even on the stack.
	for(int count = 0; count < SEGMENTS; count++) {
		seg[count].Vvb = rs.UniformRandom(0, 20);
		seg[count].Vdb = rs.UniformRandom(0, 20);
		seg[count].theta[0] = 0;
		seg[count].theta[1] = 0;
		seg[count].theta[2] = 0;
	}
	
	// Used as count variables
	int i = 0;
	int j = 0;
	
	// 1st derivative accumulator for tail segment
	// Used to calculate fitness of oscillation
	double diff[SEGMENTS] = {0};

	// Variables used to calculate phase lag fitness
	double last_local_max[SEGMENTS] = {0};
	double phase_lag_for[SEGMENTS] = {0};
	double amplitude_for[SEGMENTS] = {0};

	// good tail oscillator 0.244964 -0.476186 0.00341673 1 -0.988656 0.996945 1 0.472069 -0.00502908
	// Rename and scale evolution parameters
	double tailThreshVb = 16*0.244964;
	double tailThreshDb = 16*-0.476186;
	double tailThreshTheta = 16*0.00341673;
	double tailWeightThetaDb = 16*1;
	double tailWeightThetaVb = 16*-0.988656;
	double tailWeightDbTheta = 16*0.996945;
	double tailWeightVbTheta = 16*1;
	double tailWeightVbVb = 16*0.472069;
	double tailWeightDbDb = 16*-0.00502908;
	double bodyThreshVb = 16*t[1];
	double bodyThreshDb = 16*t[2];
	double bodyThreshTheta = 16*t[3];
	double bodyWeightThetaDb = 16*t[4];
	double bodyWeightThetaVb = 16*t[5];
	double bodyWeightDbTheta = 16*t[6];
	double bodyWeightVbTheta = 16*t[7];
	double bodyWeightVbVb = 16*t[8]*0;
	double bodyWeightDbDb = 16*t[9]*0;
	double weightVbConnection = 16*t[10];
	double weightDbConnection = 16*t[11];

	for(double trialTime = 0; trialTime <= MAXTIME; trialTime += TIMESTEP) {

		/**** SIMULATION STEP ****/
		// calculate outputs for all elements in all segments
		// from their current internal states
		double outVb[SEGMENTS];
		double outDb[SEGMENTS];
		double outTheta[SEGMENTS];

		//update tail outputs
		outVb[0] = sigmoid(seg[0].Vvb+tailThreshVb);
		outDb[0] = sigmoid(seg[0].Vdb+tailThreshDb);
		outTheta[0] = sigmoid(seg[0].theta[0]+tailThreshTheta);

		//update body outputs
		for(j = 1; j<SEGMENTS; j++){
			outVb[j] = sigmoid(seg[j].Vvb+bodyThreshVb);
			outDb[j] = sigmoid(seg[j].Vdb+bodyThreshDb);
			outTheta[j] = sigmoid(seg[j].theta[0]+bodyThreshTheta);
		}

		// affect tail accumulators
		seg[0].Vvb += (-seg[0].Vvb + tailWeightThetaVb*outTheta[0] + tailWeightVbVb*outVb[0])*TIMESTEP;
		seg[0].Vdb += (-seg[0].Vdb + tailWeightThetaDb*outTheta[0] + tailWeightDbDb*outDb[0])*TIMESTEP;
		seg[0].theta[2] = seg[0].theta[1];
		seg[0].theta[1] = seg[0].theta[0];
		seg[0].theta[0] += (tailWeightVbTheta*outVb[0] - tailWeightDbTheta*outDb[0])*TIMESTEP;

		// affect body accumulators
		for(j = 1; j<SEGMENTS; j++){
			seg[j].Vvb += (-seg[j].Vvb + bodyWeightThetaVb*outTheta[j] + bodyWeightVbVb*outVb[j] + weightVbConnection*outTheta[j-1])*TIMESTEP;	
			seg[j].Vdb += (-seg[j].Vdb + bodyWeightThetaDb*outTheta[j] + bodyWeightDbDb*outDb[j] + weightDbConnection*outTheta[j-1])*TIMESTEP;
			seg[j].theta[2] = seg[j].theta[1];
			seg[j].theta[1] = seg[j].theta[0];
			seg[j].theta[0] += (bodyWeightVbTheta*outVb[j] - bodyWeightDbTheta*outDb[j])*TIMESTEP;
		}
		

		/**** PRINT OUT DATAPOINT ****/
		// Format: trialtime angle[i] outVb[i] outDb[i]...
		if(output) {
			if(i%PRINTSTEP == 0){
				ofile << trialTime;
				for (j = 0; j<SEGMENTS; j++){
					ofile << " " << sigmoid(seg[j].theta[0]);
					ofile << " " << outVb[j];
					ofile << " " << outDb[j];
				}
				ofile << endl;
				i = 1;
			}
		}
		i++;

		/**** 1st DERIVATIVE ACCUMULATOR ****/
		for (j = 0; j < SEGMENTS; j++){
			diff[j] += fabs(sigmoid(seg[j].theta[0])-sigmoid(seg[j].theta[1]))*TIMESTEP;
		}

		/**** At T=4000, record next localmax for each segment, in sequence ****/
		if (trialTime*TIMESTEP > 4000) {
			for(int k = 0; k < SEGMENTS; k++) {
				if (last_local_max[k] == 0) {
					if (seg[k].theta[1] > seg[k].theta[0] && seg[k].theta[1] > seg[k].theta[2]) {
						last_local_max[k] = trialTime - TIMESTEP;
						amplitude_for[k] = seg[k].theta[0];
						break;
					} else {
						break;
					}
				}
			}
		}

	}
	
	//double head_oscillation = LIMIT(diff[0]/100);
	double body_oscillation = 0;
	for (j = 1; j<SEGMENTS; j++){
		body_oscillation += diff[j];
	}

	body_oscillation = LIMIT(body_oscillation/(100*(SEGMENTS-1)));
	if (body_oscillation <0.01){
		return body_oscillation;
	} else {
		double avg_phase_lag = 0;
		double phase_lag_deviation = 0;

		if (last_local_max[SEGMENTS-1] != 0) {
			// all of them oscillated, calculate phase lag phitness
			for (int k = 1; k<SEGMENTS; k++){
				phase_lag_for[k] = last_local_max[k] - last_local_max[k-1];
			}
		}

		for(int k = 1; k < SEGMENTS; k++){
			avg_phase_lag += phase_lag_for[k];
		}
		avg_phase_lag = avg_phase_lag / (SEGMENTS-1);

		for(int k = 1; k < SEGMENTS; k++){
			phase_lag_deviation += pow((phase_lag_for[k] - avg_phase_lag), 2);
		}

		return 0.01 + LIMIT(1/(phase_lag_deviation+0.000000001));
	}
}

ofstream evolution_file;

int main(int argc, char* argv[]){



if (argc > 1) {

	XG_SEED = atol(argv[1]);
	
	evolution_file.open(("evolution."+SEEDSTR+".out").c_str());
	cout.rdbuf(evolution_file.rdbuf());

	TSearch s(11);
	s.SetRandomSeed(XG_SEED);
	output = false;
	s.SetEvaluationFunction(simulate);
	s.SetSelectionMode(RANK_BASED);
	s.SetReproductionMode(GENETIC_ALGORITHM);
	s.SetPopulationSize(POPULATIONSIZE);
	s.SetMaxGenerations(GENERATIONS);
	s.SetMutationVariance(0.1);
	s.SetCrossoverProbability(0.5);
	s.SetCrossoverMode(TWO_POINT);
	s.SetMaxExpectedOffspring(1.1);
	s.SetElitistFraction(0.1);
	s.SetSearchConstraint(1);
	s.SetCheckpointInterval(5);
	s.InitializeSearch();
	s.ExecuteSearch();
	// Record best agent found to file
	ofstream best_agent;
	best_agent.open(("best_agent."+SEEDSTR+".out").c_str());
	best_agent << s.BestIndividual() << endl;
	best_agent.close();

	output = true;
	RandomState r;
	ofile.open(("best_agent_graph_hundredth."+SEEDSTR+".out").c_str());
	simulate(s.BestIndividual(), r);
	ofile.close();
	TIMESTEP = TIMESTEP/10;
	ofile.open(("best_agent_graph_thousandth."+SEEDSTR+".out").c_str());
	simulate(s.BestIndividual(), r);
	ofile.close();


	//evolution_file.close(); //This will segfault for some reason. Look into it.
		
} else {

	cout << endl;
	cout << "usage: ./command [xg_seed]" << endl;
	cout << "\t evolve and both writes to:" << endl;
	cout << "\t -best_agent.out" << endl;
	cout << "\t -best_agent_graph.out" << endl;
	cout << "\t -evolution.out" << endl;
	cout << endl;
	cout << "\t test reads from" <<endl;
	cout << "\t -best_agent.out" << endl;
	cout << "\t and writes to"<<endl;
	cout << "\t -best_agent_graph.out" << endl;
	cout << endl;

}


}
