#ifndef FROG_H
#define FROG_H

#define DIMENSION 3
#define GM 10.0   //G=1
//Mass of galaxy M = 10 solar_mass

using namespace std;
using std::setprecision;  //for decimal precision...

//========================================================================
//frog class

class frog {
private:
     ;
public:
    void leapfrog(double h, potential *Phi, star *getf);
    void drift(double h, star *posvel);
    void kick(double h, double *force, star *posvel);
 
};
//======================================================================

#endif
