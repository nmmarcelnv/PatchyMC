#ifndef _MCVOL_H
#define _MCVOL_H
#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include "Config.h"
#include "PatchyInteraction.h"
#include "rand.h"
#include "BaseParticle.h"

void MCVOL(
        BaseParticle **particles,
        int &Attempt,
        int &Acc,
        double Vmax,
        double En[2],
        double Vir[2],
        Random &rng,
        PatchyInteraction *ep
);
	

#endif  //_MCVOL_H
