/******* C Elegans Model ********
	Author: Jordan Perr-Sauer
 
 // E---
 // Warning:
 // The body starts at 0
 // J---
 // The neural network starts at 1, and goes in 2s. 
 // A lot of work to change body to start at 1.
 // J--- :)
 // E---
*********************************/


#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include <ctime>

#include <TSearch.h>
#include <random.h>
#include <VectorMatrix.h>
#include <CTRNN.h>

#include <LinearBody.h>

using namespace std;


/**** Simulation constants ****/
#define MAXTIME 300     // Simulation time
#define SEGMENTS 11     

/**** Evolution constants ****/
#define GENERATIONS 200
#define POPULATIONSIZE 20
//#define CENTER_CROSSING

/**** Convenience macro, return seed as string ****/
#define SEEDSTR ltos(XG_SEED)
string inline ltos(long x){
	char a[50];
	sprintf(a, "%ld", x);
	return a;
}


//#define LIMIT(x) (2*sigmoid(x)-1)
#define VB_i ((2*i)+1)
#define DB_i ((2*i)+2)


/**** Globals and defaults ****/
// Simulation timestep
double TIMESTEP = 0.01;

double BODYPREC = 1; // How many more times the body will integrate compared to the neural network

// TO BE IMPLEMENTED...
const double WEIGHTRANGE = 16;
const double STRETCHRECEPTORRANGE = 1;
const double BIASRANGE = 16;
// const double TIMECONSTANTRANGE 

// First argument to program. Used as seed and in file names.
long XG_SEED;
// Some internal globals
bool output = false;
bool initialization_phase = false;
ofstream physics_file;
ofstream evolution_file;
ofstream ctrnn_file;

double eval_function(TVector<double> &t, RandomState &rs)
{
	// Instantiate the neural network and body
	CTRNN net(2*SEGMENTS);      // all connections start off at 0.0
	LinearBody body(SEGMENTS);

	// Name the TVector params
    // Genotype phenotype mapping
    // Theta is the joint, angle
    
    // Body are the segments that are not the tail
	double bodyBiasVb = 16*t[1];
	double bodyBiasDb = 16*t[2];
	double bodySelfConnectionWeightVb = 16*t[3];
	double bodySelfConnectionWeightDb = 16*t[4];
	double bodyWeightThetaDb = 16*t[13];
	double bodyWeightThetaVb = 16*t[14];
	double bodyWeightDbTheta = t[15];
	double bodyWeightVbTheta = t[16];
    
    // Head is the tail
    double headBiasVb = 16*t[5];
	double headBiasDb = 16*t[6];
	double headSelfConnectionWeightVb = 16*t[7];
	double headSelfConnectionWeightDb = 16*t[8];
	double headWeightThetaDb = 16*t[9];
	double headWeightThetaVb = 16*t[10];
	double headWeightDbTheta = t[11];
	double headWeightVbTheta = t[12];

    // these two are the stretch receptors, connecting tail with last piece of the body and so on
	double weightVbConnection = 16*t[17];
	double weightDbConnection = 16*t[18];
    
    // ---- new
	double headTimeConstant = 3*(t[19]+1.2);
	double bodyTimeConstant = 3*(t[20]+1.2);

	// Configure neurons
    // 1 and 2 are the tail, V and D. 
    // 3 and 4 are the segment anterior to the tail, 
    // and so on
	net.SetNeuronBias(1, headBiasVb);
	net.SetNeuronBias(2, headBiasDb);
	net.SetNeuronTimeConstant(1, headTimeConstant);
	net.SetNeuronTimeConstant(2, headTimeConstant);
	net.SetConnectionWeight(1, 1, headSelfConnectionWeightVb);
	net.SetConnectionWeight(2, 2, headSelfConnectionWeightDb);
	for(int i = 1; i < SEGMENTS; i++) {
		net.SetNeuronBias(VB_i, bodyBiasVb);
		net.SetNeuronBias(DB_i, bodyBiasDb);
		net.SetNeuronTimeConstant(VB_i, bodyTimeConstant);
		net.SetNeuronTimeConstant(DB_i, bodyTimeConstant);
		net.SetConnectionWeight(VB_i, VB_i, bodySelfConnectionWeightVb);
		net.SetConnectionWeight(DB_i, DB_i, bodySelfConnectionWeightDb);
	}

	net.RandomizeCircuitState(0, 0);
    
#ifdef	CENTER_CROSSING         
    // stretch receptors need to be taken into account
	if(initialization_phase){
		net.SetCenterCrossing();
		t[1] = net.NeuronBias(3)/16;
		t[2] = net.NeuronBias(4)/16;
		t[5] = net.NeuronBias(1)/16;
		t[6] = net.NeuronBias(2)/16;
		return 0.0;
	}
#endif

	// Configure body
	body.nodeMass = 0.1; // KG
	body.springK = 1.5; // newtons / meter
	body.springL = 2.0 / 11; // meters
	body.springD = 1.1; // newtons / meter / second
	body.muscleK = 0.05; // newtons / radian
	body.torsionalSpringD = 0.001; // newtons / radians / second
	body.torsionalSpringK = 0.2; // newtons / radian
	body.fluidDragConstant = 100; // newtons / meter / second
	body.nodeFriction = 0.8; // newtons / meter / second
	body.Setup(TIMESTEP/BODYPREC);

	// Simulation loop
	double stretchReceptorOutput[SEGMENTS];  // angle output is restricted to [-pi, pi]
	double* COM = body.GetCOM();
	double oldx = COM[0];
	for(double trialTime = 0; trialTime <= MAXTIME; trialTime += TIMESTEP) {

		// Precompute stretch receptor outputs (between [0,1])
        // When the worm is straight, the stretch_receptor_output is 0.5
		for (int i = 0; i < SEGMENTS; i++) {
			stretchReceptorOutput[i] = body.GetBendingAngle(i)/(2*M_PI);
		}
        
		// Calculate and apply neuron external inputs from stretch receptors
		net.SetNeuronExternalInput(1, headWeightThetaVb*stretchReceptorOutput[0]);
		net.SetNeuronExternalInput(2, headWeightThetaDb*stretchReceptorOutput[0]);
		for (int i = 1; i < SEGMENTS; i++) {
			net.SetNeuronExternalInput(VB_i, bodyWeightThetaVb*stretchReceptorOutput[i] + weightVbConnection*stretchReceptorOutput[i-1]);
			net.SetNeuronExternalInput(DB_i, bodyWeightThetaDb*stretchReceptorOutput[i] + weightDbConnection*stretchReceptorOutput[i-1]);
		}

		// Update neural network with the updated inputs
		net.EulerStep(TIMESTEP);

		// Update muscle inputs using new neuron outputs
        // 0 is the Tail
        // Input to the muscle is bounded between [-1, 1]
		body.SetMuscleInput(0, (headWeightVbTheta*net.NeuronOutput(1) + headWeightDbTheta*net.NeuronOutput(2)) / 2);
		for (int i = 1; i < SEGMENTS; i++) {
			body.SetMuscleInput(i, (bodyWeightVbTheta*net.NeuronOutput(VB_i) + bodyWeightDbTheta*net.NeuronOutput(DB_i)) / 2);
		}

		// Update body with updated inputs
        // If the body needs to be integrated with more precision
		for (int i = 0; i<BODYPREC; i++){
			body.EulerStep();
		}

		if (output){
			// Write body state to physics file
			body.GetCOM();
			physics_file << COM[0] << " " << COM[1];
			nodeState* nodes = body.GetNodesPointer();
			for (int j = 0; j < body.GetNumSegments()+2; j++) {
				physics_file << "\t" << nodes[j].p.x << "\t" << nodes[j].p.y;
			}
			physics_file << endl;

			// Write neuron state to neuron file
			ctrnn_file << trialTime;
			for (int j = 1; j <= net.CircuitSize(); j++){
				ctrnn_file << "\t" << net.NeuronOutput(j);
			}
			ctrnn_file << endl;
		}

		if (body.dead)
			return 0.0;

	}
	body.GetCOM();
	if isnan(fabs(COM[0]-oldx))
		return 10000000.0;
	return (fabs(COM[0]-oldx));
	
}

void evolution_loop()
{
	// Instantiate and set up the EA
	TSearch s(20);
	s.SetRandomSeed(XG_SEED);
	output = false;
	s.SetEvaluationFunction(eval_function);
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
	RandomState r;
#ifdef CENTER_CROSSING
	initialization_phase = true;
	for (int i = 1; i<s.PopulationSize(); i++)
		eval_function(s.Individual(i), r);
#endif
	initialization_phase = false;
	s.ExecuteSearch();
    
	// Record best agent found to file
	ofstream best_agent;
	best_agent.open(("best_agent."+SEEDSTR+".out").c_str());
	best_agent << s.BestIndividual() << endl;
	best_agent.close();

    // Evaluating the best agent
	output = true;
	physics_file.open(("best_agent_hundredth_physics."+SEEDSTR+".out").c_str());
	ctrnn_file.open(("best_agent_hundredth_ctrnn."+SEEDSTR+".out").c_str());
	eval_function(s.BestIndividual(), r);
	physics_file.close();
	ctrnn_file.close();
	TIMESTEP = TIMESTEP/10;
	physics_file.open(("best_agent_thousandth_physics."+SEEDSTR+".out").c_str());
	ctrnn_file.open(("best_agent_thousandth_ctrnn."+SEEDSTR+".out").c_str());
	eval_function(s.BestIndividual(), r);
	physics_file.close();
	ctrnn_file.close();

}

// main -> Process command line arguments and invoke evolution loop
int main(int argc, char* argv[]){
  if (argc < 2){
    cout << "You need to supply an integer seed as an argument" << endl;
    return 0;
  }
	XG_SEED = atol(argv[1]);
	evolution_file.open(("evolution."+SEEDSTR+".out").c_str());
	cout.rdbuf(evolution_file.rdbuf());
	evolution_loop();
}

