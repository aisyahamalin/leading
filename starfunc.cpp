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
//Mass of galaxy M = 10 units of [10^12 M_sun]
#define PI 3.142857


//========================================================================
//definitions for star class
//========================================================================
void star::setstar(double *x, double *v) {
    for (int i = 0; i < DIMENSION; i++) {
        q[i] = x[i];   //assigning the position
        p[i] = v[i];   // assigning the velocity
    }
}//----------------------------------------------------------------------------------------------------
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
void star::printqp(ofstream& fileoutqp){
    fileoutqp << q[0] << " " << q[1] << " " << q[2] << " " << p[0] << " " << p[1] << " " << p[2] << endl; 
}
void star::printEL(ofstream& fileoutEL){
    fileoutEL << radius << " " << azi << " " << E << " " << mom[0] << " " << mom[1] << " " << mom[2] << " " << Eo << " " << tol << endl; 
}
void star::printaeTr(ofstream& fileoutaeTr){
    fileoutaeTr << a << " " << ecc << " " << b << " " << Tr << " " << x_ana << endl; 
}



//================================================================================================
//================================================================================================
//functions for returning real values
//Vc

//Mass

//


//================================================================================================
//================================================================================================


//defining the 'a' from energy: 
double star::geta(){
        double a = -GM/(2.0*E);
        return a;
}
//defining the eccentricity from 'a'
double star::geteccen(){ 
        double ecc = sqrt(1.0-(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2])/(a*GM));
        return ecc;
}
//getting the semi-minor axis b
double star::getb(){
    double b = a*sqrt(1-ecc*ecc); 
    return b;
}
//defining the radial period from 'a'
double star::getTr(){
    double Tr = 2.0*PI*sqrt(a*a*a/GM);
    return Tr; 
} 


double star::analytic(double etha[1000], double t_i){
    double M = (2*PI*t_i)/Tr; //Tr is a private member
    double RHS;
    //double x_ana;
    //double y_ana;
    for (int i=0; i<1000; i++){
        RHS = etha[i] - ecc*sin(etha[i]); //evaluating the RHS for each value of etha
        //cout << RHS << endl;

        if (abs(M-RHS)< 0.001){
            x_ana = a*(cos(etha[i])-ecc);
            //y_ana = b*sin(etha[i]);

            cout << x_ana << endl;
            return x_ana;

            break;

        }
        else {
            return 0;
            break;
        }

    }

    //return RHS;

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


//==========================================================================================
//============================================================================================
//Leapfrog method
//void function for the drift [updating the position of x, y, z each time]
void star::drift(double h) {
        for (int i = 0; i < DIMENSION; i++) {
            q[i] += (h*p[i]);
        } //position += h*velocity
}
//void function for the kick [updating the velocity using force/acceleration]
void star::kick(double h, double *force) {
    for (int i = 0; i < DIMENSION; i++) {
        p[i] += h*force[i];
    } //velocity += h*force  (force is the acceleration with point mass m)
}
//void function for the leapfrog method!!!
void star::leapfrog(double h, potential *Phi) {
    double h2 = h*0.5;
    double *force = new double[3];  //'new' operator requesting to allocate memory dynamically, here there are 3 elements...
    getforce(force, Phi);
    drift(h2);        //drift
    kick(h, force);   //kick
    drift(h2);        //drift
    
    delete [] force; //deleting memory so it won't take up space allocated by the 'new' operator
    
    E = getE(Phi); //executes the energy at each integration step
    azi = getazi();
    getcross(); //to evaluate the momentum
    radius = getr();
    tol = tolE(Phi);

    a = geta();
    ecc = geteccen();
    b = getb();
    Tr = getTr();
    //x_ana = ana(etha,t_i);
    //cout << x_ana << endl;


}
//================================================================================================
//===================================================================================================

//===============================================================================================
//==============================================================================================




// //RK4 method
// void star::runge_kutta(double h, potential *Phi){
    
//     double *force = new double[3]; //creating memory from the pointer
//     getforce(force, Phi);
//     stepA(h,force);
//     stepB(h,force);
//     stepC(h,force);
//     stepD(h,force);
    
//     //for x, y, z directions
//     for (int i = 0; i < DIMENSION; i++){
//         //position q update
//         q[i] += h/6.0 *(k1r[i] + k2r[i] + k3r[i] + k4r[i]);
//         //velocity p update
//         p[i] += h/6.0 *(k1v[i] + k2v[i] + k3v[i] + k4v[i]);
        
//         /*
//         eps[i] += L[i]*L[i] /GM ;
//         a[i] += L[i]*L[i] /(GM*(1.0-(eps[i]*eps[i]))) ;*/
//     }

//     delete [] force;

//     E = getE(Phi);
//     azi = getazi();
//     getcross(); //to evaluate the momentum
//     radius = getr();
//     tol = tolE(Phi);
//     a = geta();
//     ecc = geteccen();
//     b = getb();
//     Tr = getTr();
//     //x_ana = ana(etha,t_i);

// }
// //==================================================================================
// //in each step, loop over each array element
// //the derivative in the k1 direction
// void star::stepA(double h, double *force){
//     for (int i = 0; i < DIMENSION; i++){
//         k1r[i] = p[i];     //k1r = velocity
//         k1v[i] = force[i];     //k1v = acceleration
//     }
// }
// //the derivative in the k2 direction
// void star::stepB(double h, double *force){
//     for (int i = 0; i < DIMENSION; i++){
//         k2r[i] = p[i] + k1r[i]*h/2.0;   //velocity + k1r h/2
//         k2v[i] = force[i] + k1v[i]*h/2.0;   //accel + k1v h/2
//     }
// }
// //the derivative in the k3 direction
// void star::stepC(double h, double *force){
//     for (int i = 0; i < DIMENSION; i++){
//         k3r[i] = p[i] + k2r[i]*h/2.0;//velocity + k2r h/2
//         k3v[i] = force[i] + k2v[i]*h/2.0;  //accel + k2v h/2
//     }
// }
// //the derivative in the k4 direction
// void star::stepD(double h, double *force){
//     for (int i = 0; i < DIMENSION; i++){
//         k4r[i] = p[i] + k3r[i]*h;     //velocity + k3r h
//         k4v[i] = force[i] + k3v[i]*h;     //accel + k3v h
//     }
// }
// //================================================================================
// //=================================================================================
