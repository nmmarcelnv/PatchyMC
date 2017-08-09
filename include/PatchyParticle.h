/*******By V. Nguemaha
 *  * July 20 2016   ****/
/**
 * @brief Incapsulates a patchy particle with 2, 3, or 4 spherical patches. 
 * Used by PatchyInteraction.
 * */

#ifndef PATCHYPARTICLE_H_
#define PATCHYPARTICLE_H_

#include "BaseParticle.h"

class PatchyParticle : public BaseParticle 
{
protected:
        my_vector *_base_patches;

        void _set_base_patches();

public:
        PatchyParticle(int N_patches);
        virtual ~PatchyParticle();

        void set_positions();

};

#endif /* PATCHYPARTICLE_H_ */

