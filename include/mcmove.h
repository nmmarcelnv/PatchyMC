#ifndef _MCMOVE_H
#define _MCMOVE_H
#include "Config.h"
#include "PatchyInteraction.h"
#include "rand.h"
#include "BaseParticle.h"
#include "defs.h"

#define TRANSLATION 0
#define ROTATION 1

class MCMOVE {
friend class GCMC;

public:

MCMOVE(){

}

void step(BaseParticle **particles, PatchyInteraction *model, Random &rng,
	  int &Attempt,int &Acc, double En[2], double Vir[2], double Dmax, double Rmax=0.0, double Tprob=1.0);


private:
	void _translate_particle(BaseParticle *p, double Dmax, Random &rng);
        void _rotate_particle(BaseParticle *p, double Rmax, Random &rng);

	my_vector get_random_vector(Random &rng);
	BaseParticle **_particles_old;
	double probTranslate;
	double _delta_E;
	double Dmax;
	double Rmax;
};
#endif /* _MCMOVE_H */
