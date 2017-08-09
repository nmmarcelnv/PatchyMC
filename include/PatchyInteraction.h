#ifndef PatchyInteraction_H
#define PatchyInteraction_H
#define PI 3.141592653589793238462643f
#define PATCHY_CUTOFF 0.18f
#include "ljs.h"
#include "BaseInteraction.h"
#include "PatchyParticle.h"
#include "math.h"


class PatchyInteraction : public BaseInteraction {
    
public:

    	PatchyInteraction(BaseParticle **, Config *, double, double, int);

    	double get_total_energy(BaseParticle **particles, int NPART, int Id);

    	double compute_pair_interaction(BaseParticle *, BaseParticle *, int Id);

	int npatch(){return N_patches;}    

	double theta_max(){return cos_theta_max;}
	
	virtual ~PatchyInteraction();

private:
	// Number of patches per particle
	int N, N_patches;
   	
	//interaction energy strength between patches	
    	double patch_patch_E;

	// Patch interaction range (patch size)
	double patch_coverage, cos_theta_max;

	// Patch-patch interaction energy at the cut-off
	double patch_E_cut;	
};

#endif

