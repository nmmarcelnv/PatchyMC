#ifndef BaseInteraction_H
#define BaseInteraction_H
#include <limits>
#include "BaseParticle.h"
#include "Config.h"
#include "ljs.h"

extern double INF;
/**
 *	Abstract Base class for the BaseInteraction potential
 */
class BaseInteraction{

public:
	//Constructor
	BaseInteraction(BaseParticle **, Config *);

	//Poiter to array of all particles	
    	BaseParticle **particles;

	//Pointer to the parameter object
    	Config *conf;

	//potential energy cut off 
	double _rcut, _sqr_rcut;
	//return the total energy of the system. can be redifined by child class
    	virtual double get_total_energy(BaseParticle **particles, int NPART, int Id){};

	//pair interaction between two particles
    	virtual double compute_pair_interaction(BaseParticle *, BaseParticle *, int Id);

	//interaction energy of particle i with everyone else. 
    	double energyi(BaseParticle **particles, int NPART, int i, int jb, int Id);

	//return total energy
    	double total_energy(BaseParticle **particles, int NPART, int Id);

};
#endif   
