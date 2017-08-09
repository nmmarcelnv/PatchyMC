#include "Config.h"

Config::Config(){
}

Config::Config(int np,int ne,int nc,int ns,double T,double dr,double Rm, 
	double pt, double rho,double rc, bool shift, bool tailco):
	NPART(np), Nequil(ne),Nprod(nc),Nsamp(ns),
	T(T), DMAX(dr), RMAX(Rm), probT(pt), rho(rho), SHIFT(shift), TAILCO(tailco){

	RC = new double[2];
	RC2 = new double[2];
	ECUT= new double[2];
	BOX= new double[2];
	HBOX= new double[2];
	NPBOX = new  int[2];

	NPBOX[0] = int(np/2);
	NPBOX[1] = int(np/2);

	BOX[0] = pow( ( NPBOX[0]/rho ), 1.0/3.0 );
	BOX[1] = pow( ( NPBOX[1]/rho ), 1.0/3.0 );
	HBOX[0]=0.5*BOX[0];
	HBOX[1]=0.5*BOX[1];

	for(int ib=0; ib<2; ib++){
		RC[ib]  = std::min(25000.0, rc);
		RC2[ib] = RC[ib]*RC[ib];
	}

	BETA = 1.0/T;
	P = 0.0;
}
void Config::reset(int N1, int N2, double L1, double L2){
	NPBOX[0] = N1;
	NPBOX[1] = N2;
	BOX[0] = L1;
        BOX[1] = L2;
        HBOX[0]=0.5*BOX[0];
        HBOX[1]=0.5*BOX[1];
}

void Config::print(){
	printf(
	" Number of particles           : %d \n"
	" Temperature                   : %lf\n"
	" Pressure                      : %lf\n"
	" Box 1 Density                 : %lf\n"
	" Box 1 length                  : %lf\n"
	" Box 1 Nparticles              : %d\n"
	" Box 2 Density                 : %lf\n"
	" Box 2 length                  : %lf\n"
	" Box 2 Nparticles              : %d\n"
	" Cut off radius                : %lf\n\n"
	,NPART,T,P,(NPBOX[0]/(BOX[0]*BOX[0]*BOX[0])),BOX[0], NPBOX[0],
	(NPBOX[1]/(BOX[1]*BOX[1]*BOX[1])),BOX[1], NPBOX[1], RC[0]);
}


