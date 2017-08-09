#ifndef DEFS_H_
#define DEFS_H_

//#define PI 3.141592653589793238462643f

#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define LRACOS(x) (((x) > 1) ? (double) 0 : ((x) < -1) ? (double) PI : acos(x))

#include <stdio.h>
#include <cmath>

class my_vector;

/**
 * @brief A three-dimensional vector with built-in basic vectorial operations.
 */

class my_vector {

public:
        // Data
        double x, y, z;

        // Constructors
        my_vector(double InX, double InY, double InZ) :
                x(InX), y(InY), z(InZ) {
        }

        my_vector() :
                x(0), y(0), z(0) {
	}

        my_vector(const my_vector &p) : x(p.x), y(p.y), z(p.z) {

        }

        my_vector (double * arg) :
                x (arg[0]), y (arg[1]), z (arg[2]) {

        }

        // Operator Overloads
        inline bool operator==(const my_vector& V2) const {
                return (x == V2.x && y == V2.y && z == V2.z);
        }

        inline my_vector operator+(const my_vector& V2) const {
                return my_vector(x + V2.x, y + V2.y, z + V2.z);
        }

        inline my_vector operator-(const my_vector& V2) const {
                return my_vector(x - V2.x, y - V2.y, z - V2.z);
        }
        inline my_vector operator-() const {
                return my_vector(-x, -y, -z);
        }

        inline my_vector operator+() const {
                return my_vector(x, y, z);
        }

        inline my_vector operator/(double S) const {
                double fInv = ((double)1.) / S;
                return my_vector(x * fInv, y * fInv, z * fInv);
        }

        inline my_vector operator/(const my_vector& V2) const {
                return my_vector(x / V2.x, y / V2.y, z / V2.z);
        }

	inline double operator*(const my_vector& V2) const {
                return x * V2.x + y * V2.y + z * V2.z;
        }	

        inline my_vector operator*(double S) const {
                return my_vector(x * S, y * S, z * S);
        }
	
	inline void operator=(const my_vector& V2) {
		x = V2.x;
		y = V2.y;
		z = V2.z;
        }

        inline void operator+=(const my_vector& V2) {
                x += V2.x;
                y += V2.y;
                z += V2.z;
        }
        inline void operator-=(const my_vector& V2) {
                x -= V2.x;
                y -= V2.y;
                z -= V2.z;
        }
        inline void operator/=(double S) {
                x /= S;
                y /= S;
                z /= S;
        }

        inline void operator*=(double S) {
                x *= S;
                y *= S;
                z *= S;
        }

        inline double operator[](int i) const {
                if (i == 0)
                        return x;
                else if (i == 1)
                        return y;
                else
                        return z;
        }

        inline my_vector cross(const my_vector &V2) const {
                return my_vector(y * V2.z - z * V2.y, z * V2.x - x * V2.z, x * V2.y - y * V2.x);
        }

	inline double dot(const my_vector& V2) const {
                return x * V2.x + y * V2.y + z * V2.z;
        }    

        double norm() const {
                return x * x + y * y + z * z;
        }

        double module() const {
                return sqrt(norm());
        }
        double sqr_min_image_distance(const my_vector &V1, const double box_side) const {
                double nx = x - V1.x;
                double ny = y - V1.y;
                double nz = z - V1.z;

                nx -= rint(nx / box_side) * box_side;
                ny -= rint(ny / box_side) * box_side;
                nz -= rint(nz / box_side) * box_side;

                return nx*nx + ny*ny + nz*nz;
        }

        double sqr_distance(const my_vector &V1) const {
                return (*this - V1).norm();
        }

        double distance(const my_vector &V1) const {
                return sqrt(sqr_distance(V1));
        }

        my_vector minimum_image(const my_vector &V1, const double box_side) {
                my_vector r = *this - V1;

                r.x -= rint(r.x / box_side) * box_side;
                r.y -= rint(r.y / box_side) * box_side;
                r.z -= rint(r.z / box_side) * box_side;

                return r;
        }
        inline void normalize() {
                double fMag = (x * x + y * y + z * z);
                if (fMag == 0) return;

                double fMult = ((double)1.) / sqrtf(fMag);
                x *= fMult;
                y *= fMult;
                z *= fMult;
                return;
        }

	inline void periodic_boundaries(const double box_size){
        	if      (x < 0) x += box_size;
         	else if (x >= box_size) x -= box_size;

        	if (y < 0) y += box_size;
        	else if (y >= box_size) y -= box_size;

        	if (z < 0) z += box_size;
        	else if (z >= box_size) z -= box_size;
	}

	void print(){
	       printf( "coordinates: x = %6.3f, y = %6.3f, z = %6.3f\n", x, y, z );
	}

};

//my_vector operator*(const double S, const my_vector &v);

class my_matrix {
public:
        my_vector v1, v2, v3;

        // Constructors
        my_matrix(const my_vector nv1, const my_vector nv2, const my_vector nv3) :
                v1(nv1), v2(nv2), v3(nv3) {
        }

        my_matrix(const double n1, const double n2, const double n3, const double n4, const double n5,
                        const double n6, const double n7, const double n8, const double n9) :
                v1(my_vector(n1, n2, n3)), v2(my_vector(n4, n5, n6)), v3(my_vector(n7, n8, n9)) {
        }

        my_matrix() {

        }

        my_matrix(const my_matrix  &p) : v1(p.v1), v2(p.v2), v3(p.v3) {

        }

        // Operator Overloads
        inline bool operator==(const my_matrix &m) const {
                return (v1 == m.v1 && v2 == m.v2 && v3 == m.v3);
        }

        inline my_matrix operator+(const my_matrix& m) const {
                return my_matrix(v1 + m.v1, v2 + m.v2, v3 + m.v3);
        }

        inline my_matrix operator-(const my_matrix& m) const {
                return my_matrix(v1 - m.v1, v2 - m.v2, v3 - m.v3);
        }
        inline my_matrix operator-() const {
                return my_matrix(-v1, -v2, -v3);
        }

        inline my_matrix operator+() const {
                return my_matrix(v1, v2, v3);
        }

        inline my_matrix operator/(double S) const {
                double fInv = ((double)1.) / S;
                return my_matrix(v1 * fInv, v2 * fInv, v3 * fInv);
        }

        inline my_vector operator*(const my_vector& v) const {
                return my_vector(v1 * v, v2 * v, v3 * v);
        }

        inline my_matrix operator*(const my_matrix& m) const {
                my_matrix tm = m.get_transpose();
                return my_matrix (
                                v1 * tm.v1, v1 * tm.v2, v1 * tm.v3,
                                v2 * tm.v1, v2 * tm.v2, v2 * tm.v3,
                                v3 * tm.v1, v3 * tm.v2, v3 * tm.v3);
        }

        inline void transpone() {
                my_vector nv1 = my_vector(v1.x, v2.x, v3.x);
                my_vector nv2 = my_vector(v1.y, v2.y, v3.y);
                my_vector nv3 = my_vector(v1.z, v2.z, v3.z);
                v1 = nv1;
                v2 = nv2;
                v3 = nv3;
        }

        inline my_matrix get_transpose() const {
                return my_matrix(v1.x, v2.x, v3.x, v1.y, v2.y, v3.y, v1.z, v2.z, v3.z);
        }
        inline double determinant() const {
                return v1.x*v2.y*v3.z + v1.y*v2.z*v3.x + v1.z*v2.x*v3.y - (v3.x*v2.y*v1.z + v3.y*v2.z*v1.x + v3.z*v2.x*v1.y);
        }

        inline void orthonormalize() {
                v1.normalize ();
                v3.normalize ();
                v1 -= v3 * (v1 * v3);
                v1.normalize ();
                v2 = v3.cross (v1);
        }
};

/*
my_vector operator*(const double S, const my_vector &v) {
        return v*S;
}
*/
#endif /* DEFS_H_ */
