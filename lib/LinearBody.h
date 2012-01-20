//*********************************
//
// LinearBody simulator class header file
// A 2d vector force engine
// Jordan Perr-Sauer, July 2011
//
//
//*********************************

#pragma once

#include "VectorMatrix.h"


/*********** nodeState Struct ***********/

typedef struct {
	double x, y;
} vect;

typedef struct {
	// Physical state
	vect p; // position, Updated in Integration
	vect v; // velocity, Updated in Integration
	vect a; // acceleration, Updated in UpdateForcess
	// Force due to
	vect fd; // stoke's drag
	vect ft; // torsional spring + muscle
	vect ff; // kinetic friction
	vect fs; // structural linear springs
	// Used only for segment nodes (index 1 -> num_nodes-1)
	double muscleInput; // value between -1 and 1
	double bendingAngle; // bounded between 0 and 2PI, will not wrap around properly.
} nodeState;


// Hack to achieve modulus with doubles
double inline mod(double a, double b)
{
int result = static_cast<int>( a / b );
return a - static_cast<double>( result ) * b;
}


/*********** LinearBody Class ***********/

class LinearBody {
public:

	// Setup parameters. Must be set before simulation starts.
	double springK, springL, springD;
	double torsionalSpringK;
	double torsionalSpringD;
	double nodeFriction;
	double muscleK;
	double nodeMass;
	double fluidDragConstant;

	// Other parameters
	double timestep;
	bool dead;

	// Setup everything.
	void Setup(double timestep);

	// Constructor
	LinearBody(int segments);
	// Deconstructor
	~LinearBody();

	// Simulation Loop
	double GetBendingAngle(int segment); // in radians
	void SetMuscleInput(int segment, double input); // a value between -1 and 1
	double* GetCOM();
	int EulerStep();
	int RK4Step();
	nodeState* GetNodesPointer();
	int GetNumSegments();

protected:
	nodeState* nodes;
	int num_nodes;
	double* agent_position;

	// I split the force updates into separate functions
	// so they can be easily enabled/disabled/modified
	// by subclasses. I figure the expense of 4-5 function
	// calls per evaluation is worth the flexibility.

	// Force Helpers
	void UpdateForces();
	bool DeathCheck();
	void f_clear();
	void f_precalc();
	double* f_tmpx;
	double* f_tmpy;
	double* f_distances;
	double* f_distances_delta;
	double* f_bendingAngle_delta;

	// Forces
	void fa_linearSprings();
	void fa_torques();
	void fa_drag();
	void fa_friction();

};
