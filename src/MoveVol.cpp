
int MCVOL(
        BaseParticle **particles,
        int &Attempt,
        int &Acc,
        double Vmax,
        double En[2],
        double Vir[2],
        Random &rng,
        PatchyInteraction *ep,
        Config *conf)
{
	double volo[2], voln[2], ratio[2];

	Attempt++;
	for(int ib=0; ib<2; ib++){
		volo[ib] = pow(conf->BOX[ib],3);
		voln[ib] = volo[ib] + (2*rng()-1)*Vmax;
	}  	
  
	if (voln[0]<0 || voln[1]<0) return 1;
  
	ratio[0] = pow(voln[0]/volo[0],1/3.0);
	ratio[1] = pow(voln[1]/volo[1],1/3.0);
  	
	//double arg = -conf->BETA * conf->Pressure * (voln - vold) + N * log(vnew/vold);
  	
	if( exp(arg) <rng() ) return 1;

	double olden=calcenergy();
 
    	for (int i=0; i<conf->NPART; i++)
		idi  = particles[i]->boxId;
      		particle[i]->position.x *= ratio[idi];
		particle[i]->position.y *= ratio[idi];
		particle[i]->position.z *= ratio[idi];
	}
    
	for(int ib=0; ib<2; ib++){
		conf->BOX[ib] *= ratio[ib];
		conf->HBOX[ib] = 0.5*conf->BOX[ib];
	}
 

  	double newen=calcenergy();
  	if ( newen!=-1 && rng()<exp(conf->BETA * (olden - newen))){
  
    		Acc++;
    		energy+=newen-olden;
    		return 0;
  	}
  	else{    
      		for (int i=0; i<conf->NPART; i++){
			idi  = particles[i]->boxId;
			particle[i]->position.x /= ratio[idi];
			particle[i]->position.y /= ratio[idi];
			particle[i]->position.z /= ratio[idi];
  		}
		for(int ib=0; ib<2; ib++){ 
      			conf->BOX[ib] /= ratio[ib];
		}
    	}

    	return 2;  
}
