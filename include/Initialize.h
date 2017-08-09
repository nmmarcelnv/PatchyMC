/**
 * @file initialize.h
 * @brief Defines functions to initalize the particles
 */

#ifndef Initializer_H
#define Initializer_H

#include <math.h>
#include <stdio.h>
#include <string>
#include "rand.h"
#include <iostream>
#include "BaseParticle.h"
#include "defs.h"
#include <fstream>

#define MAX_TRIALS 100000
class Initializer{ 

private:
	double box_size;
	Random rng;
public:
	Initializer(){
		rng.setSeed(48573478);
	}
	void init_cube( BaseParticle **_particles, int n, double L);
	void init_lattice(BaseParticle **particles, int N, double rho);
	
	int lattice(BaseParticle **_particles, int N_particle, double rho, int start);
	void lattice(BaseParticle **particles, int N_particle, double rho);
	
	my_vector get_random_vector();
	void set_orientation(BaseParticle *p);
};
#endif   

