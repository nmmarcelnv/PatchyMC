#include "MCmove.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void load_restart(BaseParticle **particles, PatchyInteraction *model, Config *c, InputOutput io, int N, int nc, int npatch){
	double en[2];
	char fn[50];
	sprintf(fn, "restart.pdb");
	io.saveConfiguration(particles, c, fn, N, npatch, 1);

	printf("\nloading restart configuration! ...\n");
	io.loadConfiguration(particles, c, fn, N, npatch);
	printf("\nStarting configuration loaded! ...\n");
        for(int ib=0; ib<2; ib++){
                en[ib] = model->get_total_energy(particles, N, ib);
                printf("Loaded_Energy_Box_%d  : %6.4f\n", (ib+1), en[ib]);
        } 

	sprintf(fn, "restarti.pdb");
	io.saveConfiguration(particles, c, fn, N, npatch, 1);
        //printf("patch_coverage(cos(teta_max)): %6.4f\n",model->theta_max());
}

void MCmove::set_rand(int seed){
        iseed = seed;
        rng.setSeed(seed);
        printf("****** test random numbers ******\n");
        for (int i=1; i<=5; i++){
             printf("i,ranf() %4d %10.9f\n", i, rng() );
        }printf("\n");
}
 
MCmove::MCmove(BaseParticle **particles, PatchyInteraction *ep, InputOutput io, int ne, int nc, int nd, int nv, int ns, int fs, double dr, double dv, double rm) :
    particles(particles), model(ep), io(io), nequil(ne), nprod(nc), ndispl(nd), nvol(nv), nswap(ns), nsamp(fs), Dmax(dr), Vmax(dv), Rmax(rm)
{   
	 
    	//Allocate variables to accumulate averages for the 2 boxes    
    	current_E.resize(2);	total_E.resize(2);	eSum.resize(2);		delta_E.resize(2);
	density.resize(2);	nSamp.resize(2);	ICHP.resize(2);		CHP.resize(2);
	
	//Initialize variables to accumulate averages for the 2 boxes
    	for(int ib=0; ib<2; ib++){
		total_E[ib] = 0.0; current_E[ib] = 0.0; eSum[ib] = 0.0;	delta_E[ib] = 0.0;
		density[ib] = 0.0; nSamp[ib]     = 0;   ICHP[ib] = 0;   CHP[ib]     = 0.0;
    	}

    	probTranslate = model->conf->probT;
	NPART = model->conf->NPART;
	RMAX = Rmax;
    	_rotator = new Rotator(rng);
}

MCmove::~MCmove(){
	for(int i = 0; i <= NPART; i++) 
		delete _particles_old[i]; 
	delete [] _particles_old;

	delete _rotator;
}	

void MCmove::initialize(int restart){
	int N_patches = model->npatch();
	_particles_old = new BaseParticle *[NPART+1];

	for(int i = 0; i <= NPART; i++) {
		_particles_old[i] = new PatchyParticle(N_patches);
	}

	for(int i = 0; i <= NPART; i++) {
		//copy all particles into temporary object (to be used later)
		BaseParticle *p = _particles_old[i];
		p->index = i;
		p->copy_from(*particles[i]);
		
		//initialize positions for all patches
		if(!restart) particles[i]->set_positions();
	}

	printf("\nStarting simulation ...\n");
	for(int ib=0; ib<2; ib++){
		current_E[ib] = model->get_total_energy(particles, NPART, ib);
		nSamp[ib]++;
                printf("Initial_Energy_Box_%d  : %6.4f\n", (ib+1), current_E[ib]);
        }
	printf("patch_coverage(cos(teta_max)): %6.4f\n",model->theta_max());

	//load_restart(particles, model, model->conf, io, NPART, 100, 2);
}

void MCmove::thermalize()
{
        int ran, ncycl;
        double vir[2];
        int nmoves = ndispl + nvol + nswap;
	ncycl = nequil;
	if (ncycl != 0) std::cout <<"\nStarting equilibration...\n";

	reset();
	INIT_CHEM(0, model->conf->BETA);
	ADJUST(nAttempts, nAccepts, Dmax, Attv, Accv, Vmax, 0.5, model->conf->HBOX);

	for(int i=1; i<=ncycl; i++){
		for( int imove=1; imove<=nmoves; imove++){
			ran = rng() * (ndispl+nvol+nswap);
			if (ran < ndispl) {
				step();
			}
			else if (ran < ndispl+nvol) {
				MCVOL(particles, Attv, Accv, Vmax, &current_E[0],&vir[0], rng, model);
			}
			else{
				MCSWAP(particles, Attsw, Accsw, &ICHP[0], &CHP[0], &current_E[0],&vir[0], rng, model);
			}
		}
		if( (i % (ncycl/5) ) == 0 ) {
//			printf("======> Done %6d out of %6d %6d %6d\n", i, ncycl, model->conf->NPBOX[0], model->conf->NPBOX[1] );
			ADJUST(nAttempts, nAccepts, Dmax, Attv, Accv, Vmax, 0.5, model->conf->HBOX);
		}
	}
	printf("Equilibration complete ...\n");
}

void MCmove::production()
{
    	int ran, ib, npdbs, ncycl;
    	double vir[2], v[2];
	int nmoves = ndispl + nvol + nswap;
        ncycl = nprod;
	char filename[100];
    	sprintf(filename,"trajectories.pdb");

	npdbs = int(ncycl/100);
        if (ncycl != 0) std::cout <<"\nStaring Production ...\n";

	reset();
	INIT_CHEM(0, model->conf->BETA);    
	ADJUST(nAttempts, nAccepts, Dmax, Attv, Accv, Vmax, 0.5, model->conf->HBOX);

	for(int i=1; i<=ncycl; i++){ 
		for( int imove=1; imove<=nmoves; imove++){
			ran = rng() * (ndispl+nvol+nswap);
			if (ran < ndispl) {
				step(); 
			}
			else if (ran < ndispl+nvol) {
				MCVOL(particles, Attv, Accv, Vmax, &current_E[0],&vir[0], rng, model);
			}
			else{
				MCSWAP(particles, Attsw, Accsw, &ICHP[0], &CHP[0], &current_E[0],&vir[0], rng, model);
			}
		}

		if(i % nsamp == 0) {
			io.sample(i, &ICHP[0], &CHP[0], &current_E[0], model->conf);
			
		}

		//--save some pdbs for visualization
		if(i % npdbs == 0) {
			io.saveConfiguration(particles, model->conf, filename, NPART);
                }

		if( (i % (ncycl/5) ) == 0 ) {
//			printf("======> Done %6d out of %6d %6d %6d\n", i, ncycl, model->conf->NPBOX[0], model->conf->NPBOX[1] );
			ADJUST(nAttempts, nAccepts, Dmax, Attv, Accv, Vmax, 0.5, model->conf->HBOX);
		}
	}
	printf("Production complete ...\n\n");
	if (ncycl !=0 ){
		io.summary_moves(nAttempts, nAccepts,Dmax,Attv,Accv,Vmax,Attsw,Accsw);
		summary_energy();
		INIT_CHEM(2, model->conf->BETA);
	}
}

void MCmove::step(){
    
    	nAttempts++;
    	int i = proposeMove();
    	int ido = particles[i]->boxId;
	
    	if( Accept(ido) ){

		nAccepts++;
		current_E[ido] += delta_E[ido];
		particles[i]->position.periodic_boundaries(model->conf->BOX[ido]);

    	}else{
    
		if(move == MC_MOVE_TRANSLATION) {	
			particles[i]->position = _particles_old[i]->position;
		}else{
			particles[i]->orientation  = _particles_old[i]->orientation;
			particles[i]->set_positions();
		}
    	}
}


int MCmove::proposeMove(){
    	int ido, jb=0;
    	double Ei, Ef;
    
    	//Choose a random seed particle.
    	int pi = rng.integer(0, NPART); 
    	ido = particles[pi]->boxId;

    	//selest move type (translation or rotation?)
    	move = ( rng() < probTranslate ) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;

    	// Calculate pre-move energy.
	Ei = model->energyi(particles, NPART, pi, jb, ido);

    	if(move == MC_MOVE_TRANSLATION) {

		// Store initial state of particle and displace it.
		_particles_old[pi]->position = particles[pi]->position;
		_translate_particle(particles[pi]);

    	}else{
		//Rotation
		_particles_old[pi]->orientation  = particles[pi]->orientation;
		_rotate_particle(particles[pi]);
                particles[pi]->set_positions();

		Rmax = RMAX * ( 2.0*rng() -1.0 );
    	}

    	// Calculate post-move energy and energy change.
	Ef = model->energyi(particles, NPART, pi, jb, ido);
    	delta_E[ido] = Ef - Ei;
	//return index of moving particle
    	return pi;
}


bool MCmove::Accept(int ib)
{
    	if ( delta_E[ib] == 0  ) return true ;
    	if ( delta_E[ib] == INF) return false;
    	if ( rng() < exp(-model->conf->BETA * delta_E[ib]) ) return true;
    	else return false;
}

void MCmove::_translate_particle(BaseParticle *p) 
{ 
    	p->position.x += ( rng() -0.5 )*this->Dmax;
    	p->position.y += ( rng() -0.5 )*this->Dmax;
    	p->position.z += ( rng() -0.5 )*this->Dmax;
}

inline void MCmove::_rotate_particle(BaseParticle *p) {
        double t;
        my_vector axis;

        t = ( rng() - 0.5 ) * this->Rmax;
        axis = _rotator->get_random_vector();

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

inline void MCmove::reset(){
	nAttempts = 0;
	nAccepts  = 0;
	Attv      = 0;
	Accv      = 0;
	Attsw     = 0;
	Accsw     = 0;
}


void MCmove::summary_energy(){
    	printf("\n/****** Energy Summary ******/\n");
    	for (int ib = 0; ib <2; ib++){
                int N = model->conf->NPBOX[ib];
		if (N == 0) 
			total_E[ib] = 0;
		else
                	total_E[ib] = model->get_total_energy(particles, NPART, ib);
                double de = total_E[ib] - current_E[ib];

                if( abs( de ) > 0.000001 ) printf("**** PROBLEMS ENERGY ****\n");

                printf( "Final_Energy_Box_%d: Ent= %8.4f En= %8.4f\n"
                        ,(ib+1), total_E[ib], current_E[ib]);
   	}
}

void MCmove::INIT_CHEM(int Switch, double beta){
        if (Switch == 0) {
                for (int ib = 0; ib <2; ib++) {
                        CHP[ib]  = 0;
                        ICHP[ib] = 0;
                }
        }
        else if (Switch == 2) {
                printf("\n/****** Chemical Potential Summary ******/\n");
                for (int ib = 0; ib <2; ib++) {
                        if (ICHP[ib] !=0 )
                                CHP[ib] = -log(CHP[ib]/ICHP[ib])/beta;
                }
                printf("number of samples : %6d \n", (ICHP[0] + ICHP[1])/2);
                printf("Mu_Box_1: %6.3f \nMu_Box 2: %6.3f \n", CHP[0], CHP[1]);
        }
        else{
                printf("error chemical potential \n");
                exit(EXIT_FAILURE);
        }
}

