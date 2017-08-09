#include "mcmove.h"
void MCMOVE::MCMOVE(){

}
void MCMOVE::step(
	BaseParticle **particles, 
	int &Attempt,
        int &Acc,
	double dmax,
	double En[2],
        double Vir[2],
	Random &rng,
	PatchyInteraction *model, 
	Config *conf )
{
        double enn, eno, xn, yn, zn, viro, virn, en;
        int pi, ido, jb, move;

//      ---update the number of attempts                
        Attempt = Attempt + 1;
        jb = 0;

//      ---select a particle at random
        pi = rng.integer(0, model->conf->NPART);
        ido = model->particles[pi]->boxId;

//	---select move type (translation or rotation?)
	move = ( rng() < prob ) ? TRANSLATION : ROTATION;

//      ---calculate energy old configuration
        eno = model->energyi(model->particles, model->conf->NPART, pi, jb, ido);

	if(move == TRANSLATION){	

//      ---Store initial state of particle and displace it.
	      	particles_old[pi]->position = model->particles[pi]->position;
		_translate_particle(model->particles[pi]);

	}else{
//      ---store particle orientation
		particles_old[pi]->orientation  = model->particles[pi]->orientation;	
		_rotate_particle(model->particles[pi], Rmax, rng);
		model->particles[pi]->set_positions();

		Rmax = Rmax * ( 2.0*rng() -1.0 );
	}

//      ---calculate energy new configuration:
        enn = model->energyi(model->particles, model->conf->NPART, pi, jb, ido);

//      ---acceptance test
        double du = ( model->conf->BETA )*(enn-eno);
        if( rng() < exp(-du) || du == 0 ){
//              --accepted
                Acc = Acc + 1;
                En[ido]  += (enn - eno);  

//              ---periodic boundary condition (PBC)
		model->particles[pi]->position.periodic_boundaries(model->conf->BOX[ido]);
        }else{
		
		if(move == TRANSLATION){
//              	---restore position of particle in box
			 model->particles[pi]->position = _particles_old[pi]->position;
		}else{
//			---restore orientation of particle in box
			model->particles[pi]->orientation  = _particles_old[pi]->orientation;
			model->particles[pi]->set_positions();
        	}
	}
}

/*
bool GCMC::Accept()
{
    if ( _delta_E == 0  ) return true ;
    if ( _delta_E == INF) return false;
    if ( rng() < exp(-conf->BETA * _delta_E) ) return true;
    else return false;
}
*/
my_vector get_random_vector(Random &rng) {
        double ransq = 1.;
        double ran1, ran2;

        while(ransq >= 1) {
                ran1 = 1. - 2. * rng();
                ran2 = 1. - 2. * rng();
                ransq = ran1*ran1 + ran2*ran2;
        }

        double ranh = 2. * sqrt(1. - ransq);
        return my_vector(ran1*ranh, ran2*ranh, 1. - 2. * ransq);
}

inline void _rotate_particle(BaseParticle *p, double Rmax, Random &rng) {
        double t;
        my_vector axis;

        t = ( rng() - 0.5 ) * this->Rmax;
        axis = get_random_vector(rng);

        double sintheta = sin(t);
        double costheta = cos(t);
        double olcos = 1.0 - costheta;

        double xyo = axis.x * axis.y * olcos;
        double xzo = axis.x * axis.z * olcos;
        double yzo = axis.y * axis.z * olcos;
        double xsin = axis.x * sintheta;
        double ysin = axis.y * sintheta;
        double zsin = axis.z * sintheta;

        my_matrix R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
                                xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
                                xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

        p->orientation = p->orientation * R;
}

inline void _translate_particle(BaseParticle *p, double Dmax, Random &rng)
{
    	p->position.x += ( rng() -0.5 ) * Dmax;
    	p->position.y += ( rng() -0.5 ) * Dmax;
    	p->position.z += ( rng() -0.5 ) * Dmax;
}

