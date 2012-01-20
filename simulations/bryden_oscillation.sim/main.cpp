/******************************

Representative code for oscillatory
model in Bryden's paper.

This code evolves CTRNN parameters
to oscillate.

******************************/

// System includes
#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <ctime>

// ../EvolutionaryAgents/
#include <TSearch.h>
#include <random.h>
#include <VectorMatrix.h>

#define TRIALS 1
#define MAXTIME 1000
#define PRINTSTEP 1

#define GENERATIONS 80
#define POPULATIONSIZE 100

#define SEEDSTR ltos(XG_SEED)

using namespace std;

double TIMESTEP = 0.01;
long XG_SEED;

bool output;
ofstream ofile;

string inline ltos(long x){
	char a[50];
	sprintf(a, "%ld", x);
	return a;
}

double inline sigmoid(double x) {
	return 1/(1+exp(-x));
}

double brydenSegmentCTRNN(TVector<double> &t, RandomState &rs)
{
	// 4 weights, 2 bias, (2 gains/2 time constants?)

	// enter simulation loop
	double theta[3] = {0};
	vector<double> minmax_theta;
	vector<double> minmax_time;
	double Vvb = 0;
	double Vdb = 0;
	int i = 0;

	double diff = 0;

	for(double trialTime = 0; trialTime <= MAXTIME; trialTime += TIMESTEP) {

		// update neuron state

		double threshVb = 16*t[1];
		double threshDb = 16*t[2];
		double threshTheta = 16*t[3];
		double weightThetaDb = 16*t[4];
		double weightThetaVb = 16*t[5];
		double weightDbTheta = 16*t[6];
		double weightVbTheta = 16*t[7];
		double weightVbVb = 16*t[8];
		double weightDbDb = 16*t[9];

		double outVb = sigmoid(Vvb+threshVb);
		double outDb = sigmoid(Vdb+threshDb);
		double outTheta = sigmoid(theta[0]+threshTheta);

		Vvb += (-Vvb + weightThetaVb*outTheta + weightVbVb*outVb)*TIMESTEP;
		Vdb += (-Vdb + weightThetaDb*outTheta + weightDbDb*outDb)*TIMESTEP;
		// update theta
		theta[2] = theta[1];
		theta[1] = theta[0];
		theta[0] += (weightVbTheta*outVb - weightDbTheta*outDb)*TIMESTEP;

		double bendingAngle = sigmoid(theta[0]);

		diff += fabs(bendingAngle-sigmoid(theta[1]))*TIMESTEP;

		if(output) {
			if(i%PRINTSTEP == 0){
				ofile << trialTime<<"  "<< bendingAngle << "  "<< diff << "  " << outVb <<"  "<<outDb<<endl;
				i = 1;
			}
		}
		i++;

	}
	if (output) {double TIMESTEP = 0.1;
		ofile.close();
	}

	return diff;
}

ofstream evolution_file;

int main(int argc, char* argv[]){



if (argc > 1) {

	XG_SEED = atol(argv[1]);
	
	evolution_file.open(("evolution."+SEEDSTR+".out").c_str());
	cout.rdbuf(evolution_file.rdbuf());

		TSearch s(9);
		s.SetRandomSeed(XG_SEED);
		output = false;
		s.SetEvaluationFunction(brydenSegmentCTRNN);
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
		brydenSegmentCTRNN(s.BestIndividual(), r);
		ofile.close();
		TIMESTEP = TIMESTEP/10;
		ofile.open(("best_agent_graph_thousandth."+SEEDSTR+".out").c_str());
		brydenSegmentCTRNN(s.BestIndividual(), r);
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
