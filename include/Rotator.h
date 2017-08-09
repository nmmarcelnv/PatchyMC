/**
 * @file    Rotator.h
 * @author  Nguemaha
 */
#ifndef Rotator_H_
#define Rotator_H_
#include "rand.h"
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <string>
#include <cstdio>
#include "defs.h"

class Rotator {
private:
	Random &rng;
public:
	Rotator(Random &r);
        //
	/**
	 * @brief Generates a random vector having module 1.
	 * @return
	 */
	my_vector get_random_vector();

	/**
	 * @brief Generates a random vector inside a sphere of given radius.
	 * @param r sphere radius
	 * @return
	 */
	my_vector get_random_vector_in_sphere(double r);

	/**
	 * @brief Applies the Gram-Schmidt orthonormalization to the given matrix.
	 * @param M the matrix to be orthonormalized
	 */
	void orthonormalize_matrix(my_matrix &M);

	/**
	 * @brief Returns a matrix which generates a rotation around a random axis of a random angle, extracted between 0 and max_angle.
	 * @param max_angle
	 * @return
	 */
	my_matrix get_random_rotation_matrix(double max_angle=2*M_PI);

	/**
	 * @brief Returns a matrix which generates a rotation around a random axis of the given angle.
	 * @param angle
	 * @return
	 */
	my_matrix get_random_rotation_matrix_from_angle (double angle);

};

#endif /* UTILS_H_ */
