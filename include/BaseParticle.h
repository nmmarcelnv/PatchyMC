/*******By V. Nguemaha
 *  * July 20 2016   ****/
#ifndef _BaseParticle_H
#define _BaseParticle_H
#include "defs.h"

class BaseParticle {

public:

	BaseParticle();
        virtual ~BaseParticle();

        virtual void set_positions() {}
	virtual void copy_from(const BaseParticle &p);
	inline  void soft_copy_from(const BaseParticle *p);
	int get_index() const { return index; }
	
	/// Particle position inside the box        
    	my_vector position;

	/// BaseParticle orientational matrix
        my_matrix orientation;

	// particle index
    	int index;

	//box label in which the particle is found 0 = Box1, 1 = box 2
	int boxId;

	//number of interaction centers per particle
	int N_int_centers;

	/// Positions of all interaction centers. This array must be initialized by child classes
    	my_vector *int_centers;

};
#endif
