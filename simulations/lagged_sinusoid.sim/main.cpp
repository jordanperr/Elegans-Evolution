#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <CTRNN.h>

#include <LinearBody.h>

#define TIMESTEP 0.01
#define MAXTIME 100

ofstream ofile;

void sim(double lag) {

	LinearBody body(11);

	body.nodeMass = 0.1; // KG
	body.springK = 10; // newtons / meter
	body.springL = 2.0 / 11; // meters
	body.springD = 100; // newtons / meter / second
	body.muscleK = 0.01; // newtons / radian
	body.torsionalSpringD = 0.001; // newtons / radians / second
	body.torsionalSpringK = 0.07; // newtons / radian
	body.fluidDragConstant = 0.8; // newtons / meter / second
	body.nodeFriction = 0.5; // newtons / meter / second
	body.Setup(TIMESTEP);

	// output format
	/**

	COM[0]	COM[1]&	nodes[j].px	nodes[j].py	nodes[j+1].px...&fdp[j*2] fdp[j*2+1]	fdp[(j+1)*2] fdp[(j+1)*2+1]...
	.
	.
	.
	**/

	double *COM = body.GetCOM();
	double oldx = COM[0];
	for (double t = 0; t < MAXTIME; t += TIMESTEP){
		for (int j = 0; j<body.GetNumSegments(); j++)
			body.SetMuscleInput(j, sin(t-j*lag));
	/*			
		if (i%500 == 0){
			body.GetCOM();
			ofile << COM[0] << " " << COM[1];
			nodeState* nodes = body.GetNodesPointer();
			for (int j = 0; j < body.GetNumSegments()+2; j++) {
				ofile << "\t" << nodes[j].px << "\t" << nodes[j].py;
			}
			ofile << endl;
		}
	*/
		if(body.EulerStep()){
			break;
		}
	}
	body.GetCOM();
	ofile << lag << " " << COM[0]-oldx << endl;

}

int main(void) {
	
	ofile.open("physics.out");
	for (double lag = -2*M_PI; lag < 2*M_PI; lag += 0.01)
	{
		cout << (lag+2*M_PI)/(4*M_PI) << endl;
		sim(lag);
	}
	return 0;
}
