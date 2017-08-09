#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>

#include "Headers.h"

using namespace std;
double getL(double rho,int N){
        return (rho!=0)?pow((N/rho), 0.3333333333):1.0;
}
double getRho(double L,int N){
        return (L!=0)?double(N)/(L*L*L):0.0;
}


int main(int argc, char *argv[]){

	double T=0.0, dr=0.0, rc=1.5, L=12.5998,rho=0.05;
        int nc=1000,ne=0,fs=1,ts=1, N_particle=100, nmoves=1;
        bool tailco=0,shift=1;
        char *fname, *ofile;
        int  N_patches = 1;
	long iseed = 48573478;
        double dmax=0.09, vmax = 0.01, rmax=0.1, probT=0.5;
	int ndispl = 20, nvol = 2, nswap = 4;
	double _patchE = -1.0, patch_coverage = 1.0;
	int restart = 0; // 1 for reading restart.pdb file or 0 for initializing on lattice
	for(int i=1; i<argc; i++){
                if      (!strcmp(argv[i],"-N"))  N_particle = atoi(argv[++i]);
                else if (!strcmp(argv[i],"-rho")){
                        rho=atof(argv[++i]);
                        L=getL(rho, N_particle);
                }
                else if (!strcmp(argv[i],"-L")){
                        L=atof(argv[++i]);
                        rho=getRho(L,N_particle);
                }
                else if (!strcmp(argv[i],"-T"))   T=atof(argv[++i]);
                else if (!strcmp(argv[i],"-dr")) dmax=atof(argv[++i]);
                else if (!strcmp(argv[i],"-dR")) rmax=atof(argv[++i]);
		else if (!strcmp(argv[i],"-dv")) vmax=atof(argv[++i]);
                else if (!strcmp(argv[i],"-rc"))  rc=atof(argv[++i]);
                else if (!strcmp(argv[i],"-nc"))  nc=atoi(argv[++i]);
                else if (!strcmp(argv[i],"-ne"))  ne=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-nd"))  ndispl=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-nv"))  nvol=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-ns"))  nswap=atoi(argv[++i]);
                else if (!strcmp(argv[i],"-fs"))  fs=atoi(argv[++i]);
                else if (!strcmp(argv[i],"-ts"))  ts=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-pt"))  probT=atof(argv[++i]);
		else if (!strcmp(argv[i],"-np"))  N_patches=atoi(argv[++i]);
		else if (!strcmp(argv[i],"-pr"))  patch_coverage=atof(argv[++i]);
		else if (!strcmp(argv[i],"-res"))  restart=atoi(argv[++i]);
                else if (!strcmp(argv[i],"-h")) {
                        //usage(); exit(0);
                }else{
                        fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
                        argv[i]);
                        exit(-1);
                }
        }
	
	BaseParticle**_particles = new BaseParticle*[N_particle+1];
	for(int i = 0; i <= N_particle+1; i++) {
        //        _particles[i] = new PatchyParticle(N_patches);
        }

	// Initialise input/output object.
        InputOutput io;
	
	Config *conf = new Config(N_particle,ne,nc,fs,T,dmax,rmax,probT,rho,rc,shift,tailco);

	Initializer *Init = new Initializer;
	PatchyInteraction *ep = new PatchyInteraction(_particles, conf,_patchE, patch_coverage, N_patches);

	io.initialize(_particles, conf, N_particle, rho, N_patches, restart);
	conf->print();

	MCmove *mc_simulator = new MCmove(_particles, ep, io, ne, nc, ndispl, nvol, nswap, fs, dmax, vmax, rmax);
	mc_simulator->set_rand(iseed);
	mc_simulator->initialize(restart);
	mc_simulator->thermalize();	
	mc_simulator->production();

	conf->print();
	
	restart = 1;
	char filename[100];
        sprintf(filename,"restart.pdb");
	io.saveConfiguration(_particles, conf, filename, N_particle, N_patches, restart);	
	
	delete conf;
	delete ep;
	delete Init;
	for(int i = 0; i <= N_particle+1; i++) {
          //      delete _particles[i];
        }
	delete[] _particles;
	delete mc_simulator;
}
