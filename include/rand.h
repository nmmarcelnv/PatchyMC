/************************************************************************
*	To use, first initialize by calling RANSET( Iseed )		*
*	Then a call to RAND(Iseed) return a random number btw 0-1 	*
*									*
************************************************************************/
#ifndef RAND_H
#define RAND_H

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#define M1 714025
#define IA 1366
#define IC 150889
#define RM 1.0/M1
 

class Random{

public: 	
	double  SEED[25];
	int     *ISEED;
	double  CARRY;
	int Iseed,I24,J24;
	Random(){
		Iseed=48573478;
		RANSET();
	}

	void setSeed(int seed){
		Iseed=seed;
		RANSET();
	}

	double RANDX(){

		double  randx;
		int *p = &Iseed;
		*p = (IA*(*p)+IC) % M1 ;
		randx = *p*RM;


		if (randx < 0.0) {
			printf("*** Random number is negative *** \n");
			exit (EXIT_FAILURE);
		}

		return randx;
	}
	void RANSET(){

		int i;
		double ran;
		I24     = 24;
		J24     = 10;
		CARRY   = 0.0;
		ISEED   = &Iseed;
	/*
		get rid of initial correlations in rand by throwing     
		away the first 100 random numbers generated.
	*/
		for (i=1;i<=100;i++) ran = RANDX();

	//      initialize the 24 elements of seed
		for (i=0;i<24;i++) SEED[i] = RANDX();
	}

	double RAND(){

		double carry, rcarry, uni;
		double TWOp24=16777216.0;
		double TWOm24=1.0/TWOp24;

		uni = SEED[I24-1] - SEED[J24-1] - CARRY;

		if (uni < 0.0){
			uni = uni + 1.0;
			CARRY = TWOm24;
		}else{
			CARRY = 0.0;
		}

		SEED[I24-1] = uni;
		I24 = I24 -1;
		if (I24==0) I24 = 24;
		J24 = J24 -1;
		if (J24==0) J24 = 24;

		rcarry = uni;

		return rcarry;
	}

	
	double operator()() {
                return RAND();
        }

        //! Generate a random integer in [min, max[.
        int integer(int min, int max){
                return int(RAND() * max);
        }

};

#endif
