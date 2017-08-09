/*******By V. Nguemaha
 * July 20 2016   ****/
#include "BaseParticle.h"

BaseParticle::BaseParticle() : index(-1), boxId(-1)
{
        N_int_centers = 0;
        int_centers = NULL;
}

void BaseParticle::copy_from(const BaseParticle &p) 
{
        index = p.index;
	boxId = p.boxId;
        position = p.position;
        orientation = p.orientation;

        for(int i = 0; i < N_int_centers; i++) 
		int_centers[i] = p.int_centers[i];

}

void BaseParticle::soft_copy_from(const BaseParticle *p) 
{
	position = p->position;
        orientation = p->orientation;
}


BaseParticle::~BaseParticle() 
{
	if(int_centers != NULL) delete[] int_centers;
}
