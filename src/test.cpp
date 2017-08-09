#ifndef _Test_H
#define _Test_H

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "Headers.h"

void restart(BaseParticle **particles, PatchyInteraction *model, InputOutput io, int N, int nc, int npatch, double rho){
        double en[2];
        char fn[50];
	BaseParticle **particles = new BaseParticle *[N+1];
	for(int i = 0; i <= N+1; i++) {
                particles[i] = new PatchyParticle(npatch);
        }
	io.initialize(particles, N, rho, npatch, 0);

	for(int i = 0; i <= N; i++) {
                particles[i]->set_positions();
        }

        /*sprintf(fn, "restart.pdb");
        printf("\nloading restart configuration! ...\n");
        io.loadConfiguration(particles, fn, N, nc, npatch);
        printf("\nStarting configuration loaded! ...\n");
	*/
        for(int ib=0; ib<2; ib++){
                en[ib] = model->get_total_energy(particles, N, ib);
                printf("Loaded_Energy_Box_%d  : %6.4f\n", (ib+1), en[ib]);
        }




	for(int i = 0; i <= N+1; i++) {
                delete particles[i];
        }
        delete[] particles;    
}
#endif _Test_H 
