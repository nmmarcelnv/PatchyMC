#include "mcvol.h"
void MCVOL(
	BaseParticle **particles,
	int &Attempt,
	int &Acc,
	double Vmax,    
	double En[2],
	double Vir[2],
	Random &rng,
	PatchyInteraction *ep)

{
	double  enn[2], eno[2], f[2], volo[2], voln[2];
	double  arg, volt, dlnv, dele1, dele2, dlnv1, dlnv2, enold;
	int i, ib, idi;
//printf("%s1: En1:%6f  En2:%6f\n", __func__, En[0], En[1]);
	Attempt = Attempt + 1;

//----calculate the old energy 
        for(ib = 0; ib < 2; ib++){
                eno[ib] = ep->get_total_energy(particles, ep->conf->NPART, ib);
        }

//---calulate new volume by making random walk in ln V
	volo[0] = pow(ep->conf->BOX[0],3);
	volo[1] = pow(ep->conf->BOX[1],3);
      	volt = volo[0] + volo[1];
      	dlnv = log( volo[0]/volo[1] ) + ( rng()-0.5 ) * Vmax;
      	voln[0] = exp(dlnv)*volt / (1+exp(dlnv));
      	voln[1] = volt - voln[0];

//----rescale the box size and compute fraction change
	for(ib = 0; ib < 2; ib++){
		
        	ep->conf->BOX[ib] = pow(voln[ib], 1./3.);
         	f[ib]   = ep->conf->BOX[ib] / pow(volo[ib], 1.0/3);

		/*added following lines on Jan 13, 2017	
		ep->conf->HBOX[ib]= ep->conf->BOX[ib]/2;
                ep->conf->RC[ib]  = f[ib]*ep->conf->RC[ib];
                ep->conf->RC2[ib] = ep->conf->RC[ib]*ep->conf->RC[ib];
		*/
	}
//printf("fraction1:%6.4f, %6.4f\n", f[0], f[1]);

//---determine new coordinates
	for(i = 0; i < ep->conf->NPART; i++){
		
		idi  = particles[i]->boxId;
		particles[i]->position.x *= f[idi];
		particles[i]->position.y *= f[idi];
		particles[i]->position.z *= f[idi];
	}

//---calculate new energy
	for(int ib=0; ib<2; ib++){
		enn[ib] = ep->get_total_energy(particles, ep->conf->NPART, ib);
	}

//---acceptance test:
	dele1 = enn[0] - eno[0];
        dele2 = enn[1] - eno[1];
        dlnv1 = log(voln[0]/volo[0]);
	dlnv2 = log(voln[1]/volo[1]);
      	arg   = exp( -(ep->conf->BETA)*( dele1+dele2 -(ep->conf->NPBOX[0]+1)*dlnv1/(ep->conf->BETA) -(ep->conf->NPBOX[1]+1)*dlnv2/(ep->conf->BETA) ) );


	if (rng() < arg){
//---accepted
         	Acc = Acc + 1;
         	for(ib = 0; ib < 2; ib++){
            		En[ib]  += enn[ib]-eno[ib];
         	}
	}

	else{
//---restore the old configuration
		for(ib = 0; ib < 2; ib++){

			f[ib] = 1.0/f[ib];
            		ep->conf->BOX[ib] = ep->conf->BOX[ib]*f[ib];
            		ep->conf->HBOX[ib] = 0.5*ep->conf->BOX[ib];
 
		        /*added following lines on Jan 13, 2017		        
			ep->conf->RC[ib] = f[ib]*ep->conf->RC[ib];
                        ep->conf->RC2[ib] = ep->conf->RC[ib]*ep->conf->RC[ib];
			*/
		}
		
		for(i = 0; i < ep->conf->NPART; i++){

			idi  = particles[i]->boxId;
                	particles[i]->position.x *= f[idi];
                	particles[i]->position.y *= f[idi];
                	particles[i]->position.z *= f[idi];

		}
	}
//printf("%s2: En1:%6f  En2:%6f\n", __func__, En[0], En[1]);
}	//end of MCVOL
		

