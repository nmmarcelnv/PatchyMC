#ifndef _INPUTOUTPUT_H
#define _INPUTOUTPUT_H

#include <fstream>
#include <string>
#include <math.h>
#include "BaseParticle.h"
#include "PatchyParticle.h"
#include "Initialize.h"
#include "Config.h"

class InputOutput : public Initializer
{

private:
	//Initialization option = ( 0: initialze particle on lattice), ( 1: read restart file)
	int option; 
public:
    	// constructors.
    	InputOutput();
	InputOutput(int option_);

	//Initializer
	void initialize(BaseParticle **particles,  Config *c, int N, double rho, int npatch=0, int option=0);

    	// Load a restart configuration from a plain text file.
    	void loadConfiguration(BaseParticle **particles,  Config *c, char filename[100], int N, int npatch);	

    	// Save a restart configuration to a plain text file.
    	void saveConfiguration(BaseParticle **particles, Config *c, char filename[100], int N, int npatch=0, int restart=0);

    	// Create a VMD TcL script to set the particle view and draw a bounding box.
    	void vmdScript(const double box_size);

    	void sample(int icycl, unsigned int nMu[2], double Mu[2], double E[2], Config *conf);
    	void summary_moves(int nAtt, int nAcc, double Dmax, int Attv, int Accv, double Vmax, int Attsw, int Accsw);
   
};

#endif  /* _INPUTOUTPUT_H */
