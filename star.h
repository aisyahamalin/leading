#ifndef STAR_H
#define STAR_H

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass

using namespace std;
using std::setprecision;  //for decimal precision...

//======================================================================
//========================================================================
//specifying the star class incoporating the functions
//The object of mass M that effect the orbit of point mass m
class star {
private:
    double q[DIMENSION]; //an array with 3 elements (position)
    double p[DIMENSION]; //an array with 3 elements (velocity)
    double E;
    double k1r[DIMENSION], k2r[DIMENSION],k3r[DIMENSION], k4r[DIMENSION];
    double k1v[DIMENSION], k2v[DIMENSION],k3v[DIMENSION], k4v[DIMENSION];
    double azi;
    double mom[DIMENSION];
    double radius;
    double Eo, tol;
    double a, ecc, b, Tr;
    double RHS;
    double x_ana, y_ana;


public:
    void setstar(double *x, double *v); // a constructor
    void printcoords();                
    void getforce(double *force, potential *Phi);
    void printqp(ofstream& fileoutqp);
    void printEL(ofstream& fileoutEL);
    void printaeTr(ofstream& fileoutaeTr);

    double settime(double e);//defining the function for varying h
    double getE(potential *Phi);//function for evaluating the energy
    double getazi();
    void getcross();
    double getr();
    double setEo(double E_o);
    double tolE(potential *Phi);
    double geta();
    double getb();
    double geteccen();
    double getTr();

    //analytical function
    double analytic(double etha[1000], double t_i);

    //leapfrog function
    void leapfrog(double h, potential *Phi);
    void drift(double h);
    void kick(double h, double *force);
    
    //function for applying the RK4 method
    void runge_kutta(double h, potential *Phi);
    void stepA(double h, double *force);
    void stepB(double h, double *force);
    void stepC(double h, double *force);
    void stepD(double h, double *force);
};
//======================================================================
//=======================================================================

#endif
