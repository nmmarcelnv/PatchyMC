#include "Rotator.h"

Rotator::Rotator(Random &r) : rng(r)
{
	//empty body, but necessary!
}


my_vector Rotator::get_random_vector() {
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


my_vector Rotator::get_random_vector_in_sphere(double r) {
        double r2 = SQR(r);
        my_vector res = my_vector(r, r, r);

        while(res.norm() > r2) {
                res = my_vector(
                        2.*r*( rng() - 0.5 ),
                        2.*r*( rng() - 0.5 ),
                        2.*r*( rng() - 0.5));
        }
        return res;
}

void Rotator::orthonormalize_matrix(my_matrix &m) {
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


my_matrix Rotator::get_random_rotation_matrix_from_angle (double angle) {
        my_vector axis = get_random_vector();

        double t = angle;
        double sintheta = sin(t);
        double costheta = cos(t);
        double olcos = 1. - costheta;

        double xyo = axis.x * axis.y * olcos;
        double xzo = axis.x * axis.z * olcos;
        double yzo = axis.y * axis.z * olcos;
        double xsin = axis.x * sintheta;
        double ysin = axis.y * sintheta;
        double zsin = axis.z * sintheta;

        my_matrix R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
                                xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
                                xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

        return R;
}

my_matrix Rotator::get_random_rotation_matrix(double max_angle) {
        double t = max_angle * ( rng() - 0.5 );
        return get_random_rotation_matrix_from_angle (t);
}



