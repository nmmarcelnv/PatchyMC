#include "mcvol.h"
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

void orthonormalise(my_matrix &m) {
        double v1_norm2 = m.v1 * m.v1;
        double v2_v1 = m.v2 * m.v1;

        m.v2 -= m.v1 * (v2_v1/v1_norm2) ;

        double v3_v1 = m.v3 * m.v1;
        double v3_v2 = m.v3 * m.v2;
        double v2_norm2 = m.v2 * m.v2;

        m.v3 -=  m.v1 * (v3_v1/v1_norm2) + m.v2 * (v3_v2/v2_norm2) ;

        m.v1.normalize();
        m.v2.normalize();
        m.v3.normalize();
}

void set_orientation(BaseParticle *p, Random &rng){
    p->orientation.v1 = get_random_vector(rng);
    p->orientation.v2 = get_random_vector(rng);
    p->orientation.v3 = get_random_vector(rng);
    orthonormalise( p->orientation );
}

void MCSWAP(
	BaseParticle **particles,
	int &Attempt,
	int &Acc,
	unsigned int ICHP[2],    
	double CHP[2],
	double En[2],
	double Vir[2],
	Random &rng,
	PatchyInteraction *ep)

{
	double  enn, virn, eno, viro,  xn, yn, zn;
	double  arg, vola, vold, rhoan, rhoao, rhodn, rhodo, xo, yo, zo, dele, dtaila, dtaild;

	int o, iadd, idel, jb, idi;

	Attempt = Attempt + 1;

//===select a box at random
	if(rng()<0.5){
		iadd = 0;
		idel = 1;
	}	
	else{
		iadd = 1;
         	idel = 0;
	}

	vola = pow(ep->conf->BOX[iadd],3);
	vold = pow(ep->conf->BOX[idel],3);

//---generate a random position from a uniform normal distribution for a ghost particle
	xn = ep->conf->BOX[iadd] * rng();
      	yn = ep->conf->BOX[iadd] * rng();
      	zn = ep->conf->BOX[iadd] * rng();

//---add ghost particle to box iadd to find its energy. Note that position NPART is really never accesible by energy functions
	o = ep->conf->NPART;  
	particles[o]->position.x = xn;
        particles[o]->position.y = yn;
        particles[o]->position.z = zn;
        particles[o]->boxId = iadd;
	particles[o]->index = o; 
	
	set_orientation(particles[o], rng);
	particles[o]->set_positions();
	my_matrix orientation = particles[o]->orientation;

//---calculate energy of this particle
	jb = 0;
	enn = ep->energyi(particles, ep->conf->NPART, o, jb, iadd);

//---calculate contibution to the chemical potential:
	arg = -(ep->conf->BETA)*enn;
	CHP[iadd] = CHP[iadd] + vola*exp(arg)/double(ep->conf->NPBOX[iadd]+1);

	if( ep->conf->NPBOX[iadd] == ep->conf->NPART) CHP[iadd] = CHP[iadd] + vola*exp(arg)/double(ep->conf->NPBOX[iadd]+1);
      	ICHP[iadd] = ICHP[iadd] + 1;


//---find a random particle to delete from box b and compute its energy:
	if ( ep->conf->NPBOX[idel] == 0 ) return;
	idi = -1;  //should not be 0 or 1
      	while(idi != idel){
		o = rng.integer(0, ep->conf->NPART);    
         	idi = particles[o]->boxId;
      	}
	eno = ep->energyi(particles, ep->conf->NPART, o, jb, idel);


//---acceptance test:
      	dele = enn - eno + log( vold*(ep->conf->NPBOX[iadd]+1)/(vola*ep->conf->NPBOX[idel]) ) / (ep->conf->BETA);

	if (rng() < exp( -(ep->conf->BETA)*dele ) ){
		//---accepted
         	Acc = Acc + 1;
		ep->conf->NPBOX[iadd] = ep->conf->NPBOX[iadd] + 1;
         	particles[o]->position.x = xn;
         	particles[o]->position.y = yn;
         	particles[o]->position.z = zn;
         	particles[o]->boxId = iadd;
		particles[o]->orientation = orientation;
		particles[o]->set_positions();

         	En[iadd] = En[iadd] + enn;
         	ep->conf->NPBOX[idel] = ep->conf->NPBOX[idel] - 1;
         	En[idel] = En[idel] - eno;
	}

}	//end of MCSWAP
		

