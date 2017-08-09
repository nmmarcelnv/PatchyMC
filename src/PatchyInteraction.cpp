#include "PatchyInteraction.h"
/**
*  @param particles_: A pointer to the array of all particles
*  @param c: A pointer to a parameter object for the simulation                       
*/

PatchyInteraction::PatchyInteraction(
    BaseParticle **particles, Config *conf, double e, double r, int n) :

    BaseInteraction(particles, conf), patch_patch_E(e), patch_coverage(r), N_patches(n)

{
	patch_E_cut = 0.0;
	cos_theta_max = 1- 2.0*(patch_coverage)/N_patches;
	N = conf->NPART;
	for(int i = 0; i <= N; i++) {
                if (N_patches != 0) particles[i] = new PatchyParticle(N_patches);
        }
}


double PatchyInteraction::get_total_energy(BaseParticle **particles, int NPART, int Id){
    	if (conf->NPBOX[Id]==0) {
		//printf("In %s: empty box, returning 0 energy\n", __func__);
		return 0.0;
	}
	return total_energy(particles, NPART, Id);
}

double PatchyInteraction::compute_pair_interaction(BaseParticle *p, BaseParticle *q, int Id)
{

	my_vector r = p->position.minimum_image(q->position, conf->BOX[Id]);
	double r2 = r.norm();
        int bonds = 0;
	//square well component
	if (r2 < 1) return INF;
	if (r2 >= conf->RC2[Id]) return 0.0; 

	//The patchy interaction component
	r.normalize(); r *= 2.0;        //Added a factor 2 because the patch vectors are length 0.5.	
	for(int pi = 0; pi < p->N_int_centers; pi++) {

		my_vector ppatch = p->int_centers[pi];
		//ppatch.normalize();
	
		for(int pj = 0; pj < q->N_int_centers; pj++) {
	
			my_vector qpatch = q->int_centers[pj];	
			//qpatch.normalize();

			if( (r*ppatch >= cos_theta_max) && (-r*qpatch >= cos_theta_max)){
				bonds++;
			}
		}
	}
//printf("%s , ener = %4.2f\n", __func__, count * patch_patch_E);
	return bonds * patch_patch_E;
}

PatchyInteraction::~PatchyInteraction(){
	for(int i = 0; i <= N; i++) {
                if (N_patches != 0) delete particles[i];
        }
}
