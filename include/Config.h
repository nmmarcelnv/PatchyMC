#ifndef config_h
#define config_h 
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <map>
#define NPMAX 10000

class Config {

public:
	double 	*RC, *RC2, *ECUT, *BOX, *HBOX;
	int 	*NPBOX, NPART;
	
	int Nequil,Nprod,Nsamp; 
        double T, BETA, P, DMAX, VMAX, RMAX, rho, probT;
        bool   SHIFT, TAILCO;

	Config();
	
	Config(
		int np,		//--Total Number of Particles 
		int ne,		//--Number of Equilibration Cycles
		int nc,		//--Number of Production Cycles
		int ns, 	//--
		double T,
		double dr,
		double Rm, 
		double pt, 
		double rho,
		double rc, 
		bool shift, 
		bool tailco);
		
	void print();
	void reset(int N1, int N2, double L1, double L2);
    	virtual ~Config() {

		delete [] RC;
        	delete [] RC2;
        	delete [] ECUT;
        	delete [] HBOX;
        	delete [] BOX;
        	delete [] NPBOX;
    	}
};

#endif
