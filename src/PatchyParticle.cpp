#include <math.h>
#include <stdio.h>
#include "PatchyParticle.h"
#define HALF_ISQRT3 0.28867513459481292f

PatchyParticle::PatchyParticle(int N_patches) :
    BaseParticle()
{
   	this->N_int_centers = N_patches;
    	this->_base_patches = new my_vector[N_patches];
    	this->int_centers   = new my_vector[N_patches];

    	_set_base_patches();
}

PatchyParticle::~PatchyParticle(){
        delete [] _base_patches;
}

void PatchyParticle::_set_base_patches() 
{
        switch(this->N_int_centers) {
	case 1: {
                _base_patches[0] = my_vector(0, 1, 0);

                break;
        }
        case 2: {
                _base_patches[0] = my_vector(0,  1, 0);
                _base_patches[1] = my_vector(0, -1, 0);

                break;
        }
        case 3: {
                double cos30 = cos(M_PI / 6.);
                double sin30 = sin(M_PI / 6.);

                _base_patches[0] = my_vector( 0,      1,     0);
                _base_patches[1] = my_vector( cos30, -sin30, 0);
                _base_patches[2] = my_vector(-cos30, -sin30, 0);

                break;
        }
        case 4: {
                _base_patches[0] = my_vector(-HALF_ISQRT3, -HALF_ISQRT3,  HALF_ISQRT3);
                _base_patches[1] = my_vector( HALF_ISQRT3, -HALF_ISQRT3, -HALF_ISQRT3);
                _base_patches[2] = my_vector( HALF_ISQRT3,  HALF_ISQRT3,  HALF_ISQRT3);
                _base_patches[3] = my_vector(-HALF_ISQRT3,  HALF_ISQRT3, -HALF_ISQRT3);
                break;
        }
        default:
                printf("Unsupported number of patches %d\n", this->N_int_centers);
        }

        for(int i = 0; i < this->N_int_centers; i++) {
                _base_patches[i].normalize();
                _base_patches[i] *= 0.5;
        }
}


void PatchyParticle::set_positions() {
        for(int i = 0; i < this->N_int_centers; i++)
               this->int_centers[i] = this->orientation * _base_patches[i];
}

