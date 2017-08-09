#include "Initialize.h"
/**
 *      Initialize particles on a cubic lattice. 
 */
my_vector Initializer::get_random_vector() {
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

void orthonormalize_matrix(my_matrix &m) {
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

void Initializer::set_orientation(BaseParticle *p){
    p->orientation.v1 = get_random_vector();
    p->orientation.v2 = get_random_vector();
    p->orientation.v3 = get_random_vector();
    orthonormalize_matrix( p->orientation );
}

void Initializer::init_lattice(BaseParticle **particles, int N, double rho)
{
    int itel, n, ib;
    double del;
    /**
         place N particles on a lattice with density 'rho'
         half the number in box 1 and the other half in box 2
         n : number of particle per dimension
    */

    double L1 = pow( (N/(2*rho)), 0.3333333333 );
    del = pow(L1*L1*L1, 0.3333333333);

    n = int( pow( N, 0.3333333333 ) ) + 1;
    if (n == 0 ) n = 1;
    del = del / double(n) ;
    itel = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (ib = 0; ib < 2; ib++){
                    if (itel < N){
			
			particles[itel]->index = itel;
			particles[itel]->position.x = k*del;
                        particles[itel]->position.y = j*del;
                        particles[itel]->position.z = i*del;
                        particles[itel]->boxId = ib;
			
			particles[itel]->orientation.v1 = get_random_vector();
                    	particles[itel]->orientation.v2 = get_random_vector();
                    	particles[itel]->orientation.v3 = get_random_vector();
                     	orthonormalize_matrix( particles[itel]->orientation );

                        itel = itel + 1;
                    }
                }
            }
        }
    }
    std::cout <<"Initialisation on lattice: "<<itel << "  particles placed on a lattice"<<std::endl;
}

/// init particle positions by assigning them on a cubic grid. Modified from http://www.pages.drexel.edu/~cfa22/msim/codes/mclj_widom.c
void Initializer::init_cube(BaseParticle **_particles, int N, double L){
  int n=1;
  int ix,iy,iz;
  /* Find the lowest perfect cube, n, greater than or equal to the number of particles */
  while ((n*n*n)<N) n++;

  ix=iy=iz=0;
  for (int i=0;i<N;i++){
     BaseParticle *p = _particles[i];
     p->index = i;
     p->position.x = ((double)ix+0.5)*L/n;
     p->position.y = ((double)iy+0.5)*L/n;
     p->position.z = ((double)iz+0.5)*L/n;
     set_orientation(p);
    ix++;
    if (ix==n){
      ix=0;
      iy++;
      if (iy==n){
        iy=0;
        iz++;
      }
    }
  }
}

int Initializer::lattice(BaseParticle **_particles, int N_particle, double rho, int start) {

    int i, j, k, itel, n, count;
    double  dx, dy, dz, step, L;

    L = pow( ( N_particle/rho ), 1.0/3.0 );
    n = int( pow( N_particle, 1.0/3.0 ) ) + 1;
    step= L / double(n) ;

    if (n == 0 ) n = 1;

    itel = start;
    dx   = -step;
    count= 0;
    for (i=0; i<n; i++) {
        dx += step;
        if (dx > L) dx -= L;
        dy = -step;
        for (j=0; j<n; j++) {
            dy += step;
            if (dy > L) dy -= L;
            dz = -step;
            for (k=0; k<n; k++) {
                dz += step;
                if (dz > L) dz -= L;

                if (count < N_particle){

                    _particles[itel]->index = itel;
                    _particles[itel]->position = my_vector(dx, dy, dz);

                    _particles[itel]->orientation.v1 = get_random_vector();
                    _particles[itel]->orientation.v2 = get_random_vector();
                    _particles[itel]->orientation.v3 = get_random_vector();
                     orthonormalize_matrix( _particles[itel]->orientation );
                    count++;
                    itel += 1;
                }
            }
        }
    }
    printf("Initialisation done: %d particles placed on a lattice \n",itel-start);
    return itel-start;
}

void Initializer::lattice(BaseParticle **particles, int N_particle, double rho) {


    int N = int(N_particle/2);
    int start = 0;
    printf("Creating simulation boxes ..\n");
    for (int id = 0; id<2; id++){
        printf("starting at position .. %d\n", start);
        for (int j=start; j<N+start; j++) {
            particles[j]->boxId = id;
            //printf("id: %d\n", id);
        }
        start += lattice(particles, N, rho, start);
    }
    printf("Boxes created with a total of %d particles\n\n",start);
}

