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

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass


//========================================================================
//definitions for star class
//========================================================================
void star::setstar(double *x, double *v) {
    for (int i = 0; i < DIMENSION; i++) {
        q[i] = x[i];   //assigning the position
        p[i] = v[i];   // assigning the velocity
    }
}//----------------------------------------------------------------------------------------------------



//new getter function for q and p 
void star::getstar(double *pos, double *vel){
    for (int i = 0; i < DIMENSION; i++) {
        q[i] = pos[i];   
        p[i] = vel[i];   
    }
}



//void function for printing the coordinates to screen
void star::printcoords() {
    for (int i = 0; i < DIMENSION; i++) {
        cout << q[i] << " " << p[i] << " ";
    }   //the first column gives the position, the second one gives velocity
    cout << "\n";
}//------------------------------------------------------------------------------------------
//function for calculating the azimuthal angle
double star::getazi(){
        double azi = atan2(q[1],q[0]);
        return azi;
}//--------------------------------------------------------------------
//function for returning R
double star::getr(){
        double radius = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
        return radius;
}//------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//a function for calculating the force
void star::getforce(double *force, potential *Phi) {
    Phi->getforce(q, force);  //calling getforce with position q and force from potential class Phi
}                             // -> is the arrow operator, to reference individual members of classes
//----------------------------------------------------------------------------------------------------
//function for calculating the cross product for angular mom L
void star::getcross(){ //passing the p[i] and q[i] as arguments? no need as it's been defined private
        mom[0] = q[1]*p[2]-q[2]*p[1];
        mom[1] = q[2]*p[0]-q[0]*p[2];
        mom[2] = q[0]*p[1]-q[1]*p[0];
}//----------------------------------------------------------------------------------
//calculating the kinetic energy and potential, to get total energy--------------------------------------
double star::getE(potential *Phi) {
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        Ekin += 0.5*p[i]*p[i]; //0.5*velocity^2
    }
    return (Ekin + Phi->getpot(q));
}//----------------------------------------------------------------------------------------------------
//function for defining E_o as a private member of class star
double star::setEo(double E_o){
    return Eo = E_o; //Eo has been defined as a private in the header file
    }
//-------------------------------------------------------------------------------------------
//function for the energy tolerence
double star::tolE(potential *Phi){
    double Ekin = 0.0;
    for (int i = 0; i < 3; i++) {
        Ekin += 0.5*p[i]*p[i]; //0.5*velocity^2
    }
    double totE = Ekin + Phi->getpot(q);
    double toler = (totE - Eo)/Eo; //Eo is already a private member
    return (toler);
}//---------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
//print functions

    // E = getE(Phi); //executes the energy at each integration step
    // azi = getazi();
    // getcross(); //to evaluate the momentum
    // radius = getr();
    // tol = tolE(Phi);


void star::printqp(ofstream& fileoutqp){
    fileoutqp << q[0] << " " << q[1] << " " << q[2] << " " << p[0] << " " << p[1] << " " << p[2] << endl; 
}
void star::printEL(ofstream& fileoutEL){
    fileoutEL << radius << " " << azi << " " << E << " " << mom[0] << " " << mom[1] << " " << mom[2] << " " << Eo << " " << tol << endl; 
}
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//function for varying the time-step h, depending on r and v values of particle
double star::settime(double e){
    double mod_R = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    double mod_V = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    double dt;
    return dt = e*mod_R/mod_V;
}//----------------------------------------------------------------------------------------------------
