#include<iostream> 
#include<fstream>   
#include<string>
#include<functional>
#include<cmath>
#include<random>
#include<thread>
#include<mutex>
#include<algorithm>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "potential.h" //to import the header file where all the definitions live

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass

//=========================================================================
//definitions for potential class
//========================================================================
//--------------------------------------------------------------------------
//defining the getter functions
//a function returning the Kepler potential, depending on r
double potential::getpot(double *pos) {
        double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        return (-GM/(r));
}
//--------------------------------------------------------------------------
//a function giving the force, for each dimension x, y, z
void potential::getforce(double *pos, double *force) {
    double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    double rm3 = 1.0/(pow(r, 3.0));

    for (int i = 0; i < DIMENSION; i++) {
        force[i] = -GM*pos[i]*rm3;
    }
}