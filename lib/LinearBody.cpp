//*********************************
//
// LinearBody simulator class implementation
// Jordan Perr-Sauer, July 2011
//
// TODO:
//   - Might be nice to implement a
//     coordinate class. Would tidy
//     up the code quite a bit.
//
//   - Get muscleinput and bending
//     angle out of nodeState. That
//     was a bad idea in hindsight.
//
//*********************************

#include "LinearBody.h"
#include <stdlib.h>
#include <math.h>

// Convenient looping macros
#define ALLNODES_i int i = 0; i<this->num_nodes; i++
#define SEGMENTS_i int i = 1; i<this->num_nodes-1; i++ 



/*********** CONSTRUCTOR AND DESTRUCTOR ***********/

LinearBody::LinearBody(int segments)
{
	this->num_nodes = segments+2; // Two dummy nodes on ends.
	nodes = (nodeState *)calloc(num_nodes, sizeof(nodeState)); // Contains position, velocity, acceleration, bending angle, and muscle input
	agent_position = (double *)malloc(sizeof(double)*2); // the center of mass x,y tuple
	f_tmpx = (double *)calloc(num_nodes, sizeof(double)); // Two tmp variables for each node. Is useful
	f_tmpy = (double *)calloc(num_nodes, sizeof(double)); // in calculations to have this space
	f_distances = (double*)malloc(sizeof(double)*(num_nodes-1)); // distances between nodes (used for linear springs and drag)
	f_distances_delta = (double*)calloc(num_nodes-1, sizeof(double)); // del distances between nodes (used for dampening)
	f_bendingAngle_delta = (double*)calloc(num_nodes-1, sizeof(double)); // del bending angle (used for dampening)
	dead = false;
}

// Destructor, free everything allocated in constructor
LinearBody::~LinearBody()
{
	free(nodes);
	free(agent_position);
	free(f_tmpx);
	free(f_tmpy);
	free(f_distances);
	free(f_distances_delta);
	free(f_bendingAngle_delta);
}




/*********** ACCESSOR AND SETTER METHODS ***********/

nodeState* LinearBody::GetNodesPointer()
{
	return nodes;
}

int LinearBody::GetNumSegments()
{
	return this->num_nodes-2;
}

double* LinearBody::GetCOM()
{
	// All nodes have equivelant mass, return average of coords.
	agent_position[0] = 0.0;
	agent_position[1] = 0.0;
	for(ALLNODES_i){
		agent_position[0] += nodes[i].p.x;
		agent_position[1] += nodes[i].p.y;
	}
	agent_position[0] /= this->num_nodes;
	agent_position[1] /= this->num_nodes;
	return agent_position;
}

double LinearBody::GetBendingAngle(int segment)
{
	return nodes[segment+1].bendingAngle;
}

void LinearBody::SetMuscleInput(int segment, double input)
{
	nodes[segment+1].muscleInput = input;
}




/*********** SIMULATION API METHODS ***********/

// Setup initial positions and initialize all velocities and forces.
void LinearBody::Setup(double timestep)
{
	// Set object-wide timestep parameter.
	// Used in f_precalc() below.
	this->timestep = timestep;
	
	// Position all nodes at spring's length from origin
	for(ALLNODES_i)
		nodes[i].p.x = (-springL)*i;

	// Set the bending angle of all nodes to 180 degrees.
	// Note: This is an internal variable used to calculate
	//  torques. Do not think you can set this yourself.
	for(SEGMENTS_i)
		nodes[i].bendingAngle = M_PI;

	// Run the precalculation method to setup all other
	// parameters.
	f_precalc();

}

void LinearBody::UpdateForces()
{
	// BOOKKEEPING
	// Clean up from last force calculation and precalculate some helpful values
	f_clear();
	// Precalculate some helpful values into obect-wide variables
	f_precalc();

	// FORCE APPLICATORS
	fa_linearSprings();
	fa_torques();
	fa_drag();
	fa_friction();

	// SET ACCELERATION GIVEN NEW FORCES: a = f/m
	for (ALLNODES_i) {
		nodes[i].a.x = (nodes[i].ff.x + nodes[i].fd.x + nodes[i].ft.x + nodes[i].fs.x)/nodeMass;
		nodes[i].a.y = (nodes[i].ff.y + nodes[i].fd.y + nodes[i].ft.y + nodes[i].fs.y)/nodeMass;
	}

}

bool LinearBody::DeathCheck()
{
	if (!dead)
		for (SEGMENTS_i)
			if (fabs(f_bendingAngle_delta[i]) > M_PI)
				dead = true;
	return dead;

}

// Integrate forces into positions using a simple euler step.
// You should really keep the timestep <0.01 for this type of integration.
int LinearBody::EulerStep()
{
	// Recalculate all forces
	UpdateForces();

	// Check if we're dead. Return with error if already dead.
	if (DeathCheck())
		return 1;

	// Apply simple euler step using new forces
	for (ALLNODES_i) {
		// update velocity    a = F / m
		nodes[i].v.x += (nodes[i].a.x)*timestep;
		nodes[i].v.y += (nodes[i].a.y)*timestep;
		// update posiiton
		nodes[i].p.x += (nodes[i].v.x)*timestep;
		nodes[i].p.y += (nodes[i].v.y)*timestep;
	}

	return 0;
}

// Integrate forces into positions using RK4 step
int LinearBody::RK4Step()
{
	
	// Recalculate all forces
	double realTimestep = timestep;
	timestep /= 2;
	UpdateForces();

	// Check if dead
	if (DeathCheck())
		return 1;

	// Apply RK4 Step using new forces
	for (ALLNODES_i) {
		
	}
	
	return 0;
}



/*********** FORCE BOOKKEEPING METHODS ***********/

// Zero out forces. Used to clear memory from the previous step.
void LinearBody::f_clear()
{
	for (ALLNODES_i) {
		nodes[i].a.x = 0;  nodes[i].a.y = 0;	// acceleration
		nodes[i].fd.x = 0; nodes[i].fd.y = 0;	// drag
		nodes[i].ft.x = 0; nodes[i].ft.y = 0;	// torques
		nodes[i].ff.x = 0; nodes[i].ff.y = 0;	// friction
		nodes[i].fs.x = 0; nodes[i].fs.y = 0;	// springs
	}

}


// Precalculate distances between nodes, del distances, the bending angles of each segment, and the del bending angles.
void LinearBody::f_precalc()
{
	// precalculate f_distances between each node (headbound and tailbound)
	for (int i = 0; i<num_nodes-1; i++) {
		f_distances_delta[i] = f_distances[i];
		f_distances[i] = sqrt(pow(nodes[i].p.x - nodes[i+1].p.x, 2) + pow(nodes[i].p.y - nodes[i+1].p.y, 2));
		f_distances_delta[i] = f_distances[i] - f_distances_delta[i];
	}
	

	// precalculate angle at each segment (angle on west side of a body with head north, tail south)
	// wraparound is NOT accounted for. Be careful when designing agents.
	for (SEGMENTS_i) {

		double law_of_cosine = (
				pow(f_distances[i-1], 2)
				+pow(f_distances[i], 2)
				-(pow(nodes[i-1].p.x - nodes[i+1].p.x, 2) + pow(nodes[i-1].p.y - nodes[i+1].p.y, 2))
			) / (2*f_distances[i-1]*f_distances[i]);

		//cout << law_of_cosine << endl;
		double angle;
		if (law_of_cosine <= -1){
			angle = M_PI;
		} else if (law_of_cosine >= 1){
			angle = 0;
		} else {
			angle = acos(law_of_cosine);
		}

		// use cross product to determine which side of the body we're looking for
		// http://stackoverflow.com/questions/3461453/determine-which-side-of-a-line-a-point-lies
		bool leftSide = ((nodes[i].p.x - nodes[i-1].p.x)*(nodes[i+1].p.y-nodes[i-1].p.y) - (nodes[i].p.y - nodes[i-1].p.y)*(nodes[i+1].p.x-nodes[i-1].p.x)) > 0;
		if (leftSide) {
			angle = 2*M_PI-angle;
		}

		// finally, set the bending angle
		f_bendingAngle_delta[i-1] = nodes[i].bendingAngle;
		nodes[i].bendingAngle = angle;
		f_bendingAngle_delta[i-1] = nodes[i].bendingAngle - f_bendingAngle_delta[i-1];
	}
}





/*********** FORCE APPLICATION METHODS ***********/

void LinearBody::fa_linearSprings()
{
	// Apply force of linear springs in two passes
	// headbound
	for (int i = 0; i<num_nodes-1; i++) {
		// calculate magnitude of spring force
		double springForce = springK*(springL - f_distances[i]);
		double dampForce = springD*(f_distances_delta[i]/timestep);
		//cout << timestep << endl;
		// calculate unit vector going from i+1 to i
		f_tmpx[i] = (nodes[i].p.x - nodes[i+1].p.x) / f_distances[i];
		f_tmpy[i] = (nodes[i].p.y - nodes[i+1].p.y) / f_distances[i];
		// multiply spring force by the unit vector
		f_tmpx[i] *= (springForce-dampForce);
		f_tmpy[i] *= (springForce-dampForce);
		// apply leftbound spring forces to force array
		nodes[i].fs.x += f_tmpx[i];
		nodes[i].fs.y += f_tmpy[i];
	}
	// tailbound
	for (int i = 1; i<num_nodes; i++) {
		// apply negative of force in tmp array (index i-1) to point i
		nodes[i].fs.x -= f_tmpx[i-1];
		nodes[i].fs.y -= f_tmpy[i-1];
	}
}



void LinearBody::fa_torques()
{
	double hx, tx, hy, ty;
	// Apply Torsional spring and muscle torques
	for (SEGMENTS_i) {

		// hooke's law
		double springTorque = ((nodes[i].bendingAngle-M_PI)/(M_PI))*torsionalSpringK;
		// velocity dampening
		double springDampTorque = (f_bendingAngle_delta[i-1]/timestep)*torsionalSpringD;
		// muscle definition
		double muscleTorque = nodes[i].muscleInput*muscleK;
		// sum them together into one common torque
		double torque = springTorque - springDampTorque + muscleTorque;

		// Calculate leftward normal unit vectors for torque at adjacent nodes
		// headbound
		hx = -(nodes[i-1].p.y-nodes[i].p.y) / f_distances[i-1];
		hy =  (nodes[i-1].p.x-nodes[i].p.x) / f_distances[i-1];
		// tailbound
		tx = -(nodes[i].p.y-nodes[i+1].p.y) / f_distances[i];
		ty =  (nodes[i].p.x-nodes[i+1].p.x) / f_distances[i];

		// Apply torque forces to adjacent nodes
		double headTorque = torque/f_distances[i-1];
		double tailTorque = torque/f_distances[i];
		nodes[i-1].ft.x += (headTorque)*hx;
		nodes[i-1].ft.y += (headTorque)*hy;
		nodes[i+1].ft.x += (tailTorque)*hx;
		nodes[i+1].ft.y += (tailTorque)*hy;

		// Apply torque reactive force to segment node
		nodes[i].ft.x -= (headTorque)*hx + (tailTorque)*hx;
		nodes[i].ft.y -= (headTorque)*hy + (tailTorque)*hy;
	}
}


void LinearBody::fa_drag()
{
	double vx, vy, nx, ny, proj;
	// Apply segment drag
	// Henceforth "bar" := the line between node[i] and node[i+1]
	for (int i = 0; i < num_nodes-1; i++) {
		// find normal to bar facing left using cross product, keeping the length (you'd multiply it back in anyway)
		nx = -(nodes[i].p.y-nodes[i+1].p.y);
		ny =  (nodes[i].p.x-nodes[i+1].p.x);
		// find velocity vector of midpoint of bar
		vx = (nodes[i].v.x + nodes[i+1].v.x)/2;
		vy = (nodes[i].v.y + nodes[i+1].v.y)/2;
		// find projection of velocity vector on normal (dot product)
		proj = nx*vx + ny*vy;
		// divide this magnitude in half, as each node feels 1/2 of the drag
		proj /= 2;
		proj *= fluidDragConstant;
		// apply half of drag force to each adjacent node
		nodes[i].fd.x -= nx*proj;
		nodes[i].fd.y -= ny*proj;
		nodes[i+1].fd.x -= nx*proj;
		nodes[i+1].fd.y -= ny*proj;
	}

}

void LinearBody::fa_friction()
{
	// Apply kinetic friction to all nodes
	for (ALLNODES_i) {
		nodes[i].ff.x = -nodeFriction*nodes[i].v.x;
		nodes[i].ff.y = -nodeFriction*nodes[i].v.y;
	}
}

