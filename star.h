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
    double L[DIMENSION]; //an array for the angular momentum
    double k1r[DIMENSION], k2r[DIMENSION],k3r[DIMENSION], k4r[DIMENSION];
    double k1v[DIMENSION], k2v[DIMENSION],k3v[DIMENSION], k4v[DIMENSION];
    double azi;
    double mom[DIMENSION];
    double radius;
    double Eo;
    double tol;

    //double eps[DIMENSION];
    //double a[DIMENSION];


public:
    void setstar(double *x, double *v); // a constructor
//a getter for frog
    void getstar(double *pos, double *vel);


    void printcoords();                
    void getforce(double *force, potential *Phi);
    void printqp(ofstream& fileoutqp);
    void printEL(ofstream& fileoutEL);

    double settime(double e);//defining the function for varying h
    double getE(potential *Phi);//function for evaluating the energy
    double getazi();
    void getcross();
    double getr();
    double setEo(double E_o);
    double tolE(potential *Phi);
};
//======================================================================
//=======================================================================

#endif
