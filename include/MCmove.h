#ifndef MCmove_H
#define MCmove_H
#include <math.h>
#include "rand.h"
#include "BaseParticle.h"
#include "PatchyParticle.h"
#include "PatchyInteraction.h"
#include "mcvol.h"
#include "mcswap.h"
#include "adjust.h"
#include "InputOutput.h"
#include "defs.h"
#include "Rotator.h"
#include <vector>

#define MC_MOVE_TYPES 2
#define MC_MOVE_TRANSLATION 0
#define MC_MOVE_ROTATION 1

class MCmove {

public:
    	PatchyInteraction *model;
	BaseParticle **particles;
    	Rotator *_rotator;
	InputOutput io;
    	Random rng;

    	MCmove(BaseParticle **particles, PatchyInteraction *ep, InputOutput io, int ne, int nc, int nd, int nv, int ns, int fs, double dr, double dv, double rmax);
    	~MCmove();
    	void step();
   	bool Accept(int ib);
    	void thermalize();
    	void production();
    	void summary_energy();
    	void initialize(int restart=0);
    	void reset();	
    	void set_rand(int seed);
    	void INIT_CHEM(int Switch, double beta);	
    
	double Dmax, Vmax, Rmax, RMAX;
    
	int NPART;
    int iseed; 
    int nAttempts;
    int nAccepts;
   
    int nvol;
    int nswap;
    int ndispl, nequil, nprod;
    int Attv, Accv, Attsw, Accsw;        

    int nRotations; 
    int move, nsamp;
     
    
    BaseParticle **_particles_old;
    double probTranslate;
    double _delta_E, _current_E, esum;
    std::vector<double> current_E, total_E, density, delta_E, eSum, CHP;
    std::vector<unsigned int> nSamp, ICHP;

    void _rotate_particle(BaseParticle *p);
    void _translate_particle(BaseParticle *p);

    int proposeMove();
   
};
#endif
