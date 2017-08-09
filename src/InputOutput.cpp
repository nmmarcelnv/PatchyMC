#include "InputOutput.h"

InputOutput::InputOutput() {}
InputOutput::InputOutput(int option_) : option(option_) {}

void InputOutput::initialize(BaseParticle **particles,  Config *c, int N, double rho, int npatch, int option) {
	char filename[50];
	sprintf(filename, "restart.pdb");
	switch (option){
		case 0: {
			lattice(particles, N, rho); 
			break;
		}
		case 1: { 
			loadConfiguration(particles, c, filename, N, npatch);
			break;
		}
		default:
			printf("Error: %s, Invalid initialization option\n", __func__);

	}
}


void InputOutput::saveConfiguration(BaseParticle **particles, Config *c, char filename[100], int N, int npatch, int restart){
	FILE *fp;	
	if (restart) {
		fp = fopen(filename, "w");
		fprintf(fp, "%d %d\n", N, npatch);
		fprintf(fp, "%f %f\n", c->BOX[0], c->BOX[1] );
	}else{
		fp = fopen(filename, "a");
	}
        int i; double s=1.0;
        for (i=0; i<N; i++){
		BaseParticle *p = particles[i];
                fprintf(fp,"%4s%7d%3s%6s%6d    %8.4lf%8.4lf%8.4lf\n",
                       "ATOM", i+1, "N", "MET", p->boxId, p->position.x, p->position.y, p->position.z);
		int j=0;
		while(j < p->N_int_centers){
			fprintf(fp,"%4s%7d%3s%6s%6d    %8.4lf%8.4lf%8.4lf\n",
				"ATOM", i+1, "H", "MET", p->boxId, 
				(p->position.x + s*p->int_centers[j].x), 
				(p->position.y + s*p->int_centers[j].y),
				(p->position.z + s*p->int_centers[j].z)); 
			j++;
		}
        }
	fclose (fp);
}

void InputOutput::loadConfiguration(BaseParticle **particles, Config *c, char filename[100], int N, int npatch){
	std::ifstream inFile;
	std::string tmp;
	int n, id, npart, Npatch, nbox1=0, nbox2=0;
	double L1, L2;

	inFile.open(filename);
	if (!inFile.good()){printf("Error in %s: Invalid restart file!\n", __func__); exit(EXIT_FAILURE);} 
	
	//save time step, total number of particle and number of patch per particle, 
	inFile >> npart >> Npatch;
	if (N != npart or Npatch != npatch){ 
		printf("Error in %s: num particle or number patch don't match!\n", __func__); exit(EXIT_FAILURE);
	}

	//save box lengths to temp variable
	inFile >>L1 >> L2;

	for (int i=0; i<N; i++){
		BaseParticle *p = particles[i];
		
		//throw unwanted fields
		inFile >> tmp >> n >> tmp >> tmp
		
		//load particle boxId and coordinates
		>> p->boxId >> p->position.x >> p->position.y >> p->position.z; 
		if (p->boxId == 0) nbox1++; 
		if (p->boxId == 1) nbox2++;
		//load patch interaction centers relative to particle center.
		for(int j=0; j<npatch; j++){
			double px, py, pz;
			inFile >> tmp >> n >> tmp >> tmp>> id >> px >> py >> pz;
			p->int_centers[j].x = px - p->position.x;
			p->int_centers[j].y = py - p->position.y;
			p->int_centers[j].z = pz - p->position.z;
		}
	}
	if(nbox1+nbox2 != N){
		printf("Error in %s: num particle inconsistency!\n", __func__); 
		exit(EXIT_FAILURE);
	}else{ 
		//reset the configuration parameter
		c->reset(nbox1, nbox2, L1, L2);
		printf("successfully loaded configuration:\n" 
			"Placed %d particles in Box1, and %d particles in Box2 \n", nbox1, nbox2);
	}
	
	inFile.close();
}


void InputOutput::sample(int icycl, unsigned int nMu[2], double Mu[2], double E[2], Config *conf){
	double v[2], rho[2], mu[2];
	char filename[100];
	sprintf(filename,"data-N%d-T%3.2f", conf->NPART, conf->T);
        FILE * pFile;
        pFile = fopen(filename, "a");

       	v[0]   = pow(conf->BOX[0], 3);	v[1]   = pow(conf->BOX[1], 3);
        rho[0] = conf->NPBOX[0]/v[0] ;	rho[1] = conf->NPBOX[1]/v[1] ;

	for(int ib = 0; ib <2; ib++){
		if(nMu[ib] != 0) mu[ib]  = -log(Mu[ib]/nMu[ib])/conf->BETA;
	}

        fprintf (pFile,"%10d %8d %8.4f %8.4f %8.4f %8.4f %8d %8.4f %8.4f %8.4f %8.4f\n",
		 icycl, conf->NPBOX[0], v[0], rho[0], E[0], mu[0], conf->NPBOX[1], v[1], rho[1], E[1], mu[1]);
        fclose (pFile);
}

void InputOutput::summary_moves(int nAtt, int nAcc, double Dmax, int Attv, int Accv, double Vmax, int Attsw, int Accsw){
	printf("/****** Move Summary ******/\n");
	if(nAtt !=0 ){
                printf(
                "Number of att. to displ. a part. : %d\n"
		"Success                          : %d (= %5.4f)\n"
                "Maximum displacement adjusted to : %5.4f\n",
                nAtt, nAcc, 1.0*nAcc/nAtt, Dmax );
        }
	if(Attv !=0 ){
                printf("\n"
                "Number of att. to chan. vol.     : %d\n"
		"Success                          : %d (= %5.4f)\n"
                "Maximum volume adjusted to       : %5.4f\n",
                Attv, Accv, 1.0*Accv/Attv,Vmax );
        }
        if(Attsw !=0 ){
                printf("\n"
                "Number of att. to exchan. part.  : %d\n"
		"Success                          : %d (= %5.4f)\n"
		"Particle exchange not adjustable :\n",
                Attsw, Accsw, 1.0*Accsw/Attsw);
        }
}

void InputOutput::vmdScript(const double boxSize)
{
    FILE *pFile;

    pFile = fopen("vmd.tcl", "w");

    // Turn on lights 0 and 1.
    fprintf(pFile, "light 0 on\n");
    fprintf(pFile, "light 1 on\n");
    fprintf(pFile, "light 2 off\n");
    fprintf(pFile, "light 3 off\n");

    // Position the stage and axes.
    fprintf(pFile, "axes location off\n");
    fprintf(pFile, "stage location off\n");

    // Set orthographic projection.
    fprintf(pFile, "display projection orthographic\n");

    // Set drawing method to van der Waals radius.
    fprintf(pFile, "mol modstyle 0 0 VDW 1 30\n");

    // Set sensible atom radius.
    fprintf(pFile, "set sel [atomselect top \"name X\"]\n");
    fprintf(pFile, "atomselect0 set radius 0.5\n");

    // Set default particle to blue.
    fprintf(pFile, "color Name X blue\n");

    // Turn off depth cue.
    fprintf(pFile, "display depthcue off\n");

    // Define box boundaries.
    fprintf(pFile, "set minx 0\n");
    fprintf(pFile, "set maxx %5.4f\n", boxSize);
    fprintf(pFile, "set miny 0\n");
    fprintf(pFile, "set maxy %5.4f\n", boxSize);
    fprintf(pFile, "set minz 0\n");
    fprintf(pFile, "set maxz %5.4f\n", boxSize);

    // Set colours.
    fprintf(pFile, "draw materials off\n");
    fprintf(pFile, "draw color white\n");

    // Draw cube edges.
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$maxx $miny $minz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$minx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $minz\" \"$minx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\"\n");
    fprintf(pFile, "draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\"\n");
    fprintf(pFile, "draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\"\n");

    // Rotate box.
    fprintf(pFile, "rotate x by -60\n");
    fprintf(pFile, "rotate y by -30\n");
    fprintf(pFile, "rotate z by -15\n");


    fclose(pFile);
}


