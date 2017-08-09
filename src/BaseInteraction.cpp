#include "BaseInteraction.h"
#include <stdio.h>
double INF = std::numeric_limits<double>::infinity();
BaseInteraction::BaseInteraction(BaseParticle **p, Config *c) :
    	particles(p), conf(c), _rcut(c->RC[0])
{   	
	
}


double BaseInteraction::total_energy(BaseParticle **particles, int NPART, int Id)
{
    int i, jb;
    double tot_energi = 0.0;
    
    for (i = 0; i < NPART-1; i++){
	if ( particles[i]->boxId == Id ){
        	jb = i + 1;
        	tot_energi += energyi(particles, NPART, i, jb, Id);
	}
    }
    return tot_energi;
}


double BaseInteraction::energyi(BaseParticle **particles, int NPART, int i, int jb, int Id)
{
    double run_energi = 0.0;

    for (int j=jb; j<NPART; j++){
	  if ( particles[j]->boxId == Id ){
        	if (j != i){
            		run_energi += compute_pair_interaction(particles[i], particles[j], Id);
		}
          }
    }
    return run_energi;
}


double BaseInteraction::compute_pair_interaction(BaseParticle *p1, BaseParticle *p2, int Id)
{
    printf("Error in BaseInteraction: Virtual function BaseInteraction::compute_pair_energy() must be defined.\n");
}

