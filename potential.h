#ifndef POTENTIAL_H
#define POTENTIAL_H

//=========================================================================
//specifying the potential class incorporating both getpot and getforce
//this will give the potential and force
class potential {
private:
    ;
public:
    double getpot(double *pos); //a getter function
    void getforce(double *pos, double *force); //a getter function
};
//==========================================================================

#endif