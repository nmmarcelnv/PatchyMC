#include "defs.h"
#include <stdio.h>

using namespace std;

//simple testing of method in my_vector class
int main(int argc, char * argv[]){

	my_vector v;
	v.print();

	v = my_vector(1.0, 2.0, -1.0);
	v.print();

	my_vector u = v;
	u.print();

	my_vector w = u + v;
	w.print();

	my_vector minus_w = -w;
	minus_w.print();

	my_vector a = w/2.0;
	a.print();

	my_vector ratio = w/w;
	ratio.print();

	my_vector b = my_vector(0, 1.0, -1.0);
	my_vector cross = v.cross(b);
	cross.print();

	double Nv = v.norm();
	double mv = v.module();
	printf("Norm v: %4.2f %4.2f \n", Nv, mv );

	v.normalize();
	v.print();
}
