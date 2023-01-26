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
#include "star.h" //to import the header file where all the definitions live
#include "frog.h"

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass

//==========================================================================================
//============================================================================================
//Leapfrog method
//void function for the drift [updating the position of x, y, z each time]
void frog::drift(double h, star *posvel) {
    posvel->getstar(pos, vel);
        for (int i = 0; i < DIMENSION; i++) {
            q[i] += (h*p[i]);
        } //position += h*velocity
}
//void function for the kick [updating the velocity using force/acceleration]
void frog::kick(double h, double *force, star *posvel) {
    posvel->getstar(pos, vel);
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += h*force[i];
    } //velocity += h*force  (force is the acceleration with point mass m)
}

//void function for the leapfrog method!!!
void frog::leapfrog(double h, potential *Phi, star *getf) {
    double h2 = h*0.5;
    double *force = new double[3];  //'new' operator requesting to allocate memory dynamically, here there are 3 elements...

    getf->getforce(force, Phi);

    drift(h2);        //drift
    kick(h, force);   //kick
    drift(h2);        //drift
    
    delete [] force; //deleting memory so it won't take up space allocated by the 'new' operator
    




    
    // E = getE(Phi); //executes the energy at each integration step
    // azi = getazi();
    // getcross(); //to evaluate the momentum
    // radius = getr();
    // tol = tolE(Phi);

}
//================================================================================================
//===================================================================================================

