#ifndef _MCSWAP_H
#define _MCSWAP_H
#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "rand.h"
#include "Config.h"
#include "defs.h"
#include "PatchyInteraction.h"
#include "BaseParticle.h"

void set_orientation(BaseParticle *p, Random &rng);
my_vector get_random_vector(Random &rng);
void orthonormalise(my_matrix &m);

void MCSWAP(
        BaseParticle **particles,
        int &Attempt,
        int &Acc,
        unsigned int ICHP[2],
        double CHP[2],
        double En[2],
        double Vir[2],
        Random &rng,
        PatchyInteraction *ep
);	

#endif  //_MCSWAP_H
