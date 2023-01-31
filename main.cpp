#include<iostream>  //for cout, cin
#include<fstream>   //for writing and reding from files
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


using namespace std;
using std::setprecision;  //for decimal precision...

//THE MAIN FUNCTION
int main() {
    double *x = new double[DIMENSION];  //value of pointer allocated dynamical memory
    double *v = new double[DIMENSION];  //to create a memory (of actual values) in the heap
    //the values of each array
    x[0] = 9.0; x[1] = 0.0; x[2] = 0.0; //note: our sun is about 8kpc from centre
    v[0] = 0.0; v[1] = 1.0; v[2] = 0.0; //velocities in units of ~240 km/s

    star pedro; //calling it pedro //intantiating the star class with pedro
    pedro.setstar(x, v); //setting the positions x and velocities v for pedro //instantiating the object
    
    pedro.printcoords(); //printing the initial values to terminal



    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //defining the step-length //which is a constant throughout code
    double h = 1.0e-6;  //unit of timeperiod is: 1 Gyr //here, time-step h is 1 000 years
    //note: my original 1.0e-6 == 1000 years
    // 1.0e-5 == 10 000 years
    //-----------------------------------------------------------------
    //defining the parameter constant e
    double e = 1.0e-7; //1.0e-7; //my original 1.0e-7
    double dt;
    //-----------------------------------------------------------------
    //-----------------------------------------------------------------
    //defining the original energy
    double kin = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); 
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    double pot = -GM/r;
    double E_o = kin + pot;
    double Eo;
    //-----------------------------------------------------------------

    cout << E_o << " " << r << endl;

    potential Phi; //calling it Phi //making an instance of potential called phi


//------------------------------------------------------------------------------------
//defining the eta array [0:2pi]
    double etha[1000]; //one thousand elements in array
    etha[0] = 0; //setting the zeroth element
    for (int i=0; i<999; i++){
        etha[i+1] = etha[i] + 0.001*2*PI;
        //cout << etha[i] << endl; 
    }
//---------------------------------------------------------------------------------------






    //LEAPFROG=========================================================================================
    ofstream myfile ("coor_frog.dat", std::ios_base::trunc);
    ofstream file ("EL_frog.dat", std::ios_base::trunc);
    ofstream anotherfile ("aeTr_frog.dat", std::ios_base::trunc);

    ofstream timefile("time_frog.dat", std::ios_base::trunc);

    for (int i = 0; i < 1000; i++) {   //one thousand times //how many times you print to file
        for (int ii = 0; ii < 100000; ii++) { //for each one time out of a thousand, do it 100,000 times
                                             //actually implementing the leapfrog //if using h=1 000 yr --> means 100,000,000 years
            dt = pedro.settime(e);
            //declaring the E_o inside integrator
            Eo = pedro.setEo(E_o);
            //the leapfrog for that particular star
            pedro.leapfrog(h, &Phi);
            }
 
        //pedro.printcoords();         //printing the new coordinates after applying the leapfrog
        pedro.analytic(etha,h*i*100000); 

        pedro.printqp(myfile);
        pedro.printEL(file);
        pedro.printaeTr(anotherfile);
        timefile << i*100000*h << endl;



        }

        
    myfile.close();
    file.close();
    anotherfile.close();
    timefile.close();
    //=======================================================================================
    // //RUNGE-KUTTA========================================================================================
    // pedro.setstar(x, v); //setting the positions x and velocities v for pedro //instantiating the object again for Runge-Kutta
    // delete [] x; //deleting the allocated memory
    // delete [] v; //deleting from the heap
    // //-----------------------------------------------------------------------------
    // ofstream myfile2 ("coor_RK.dat", std::ios_base::trunc);
    // ofstream file2 ("EL_RK.dat", std::ios_base::trunc);
    // ofstream anotherfile2 ("aeTr_RK.dat", std::ios_base::trunc);

    // for (int i = 0; i < 1000; i++) {   //one thousand times //how many times you print to file //sampling 1000 times
    //     for (int ii = 0; ii < 100000; ii++) { //for each one time out of a thousand, do it 100,000 times
    //                                          //actually implementing the integrator
    //         dt = pedro.settime(e);
    //         //declaring the E_o inside integrator
    //         Eo = pedro.setEo(E_o);
    //         //the actual integrator
    //         pedro.runge_kutta(h, &Phi);
    //         }

    //     //pedro.printcoords();         //printing the new coordinates after applying the leapfrog

    //     //pedro.ana(etha,h*i*100000); 



    //     pedro.printqp(myfile2);
    //     pedro.printEL(file2);
    //     pedro.printaeTr(anotherfile2);

    //     }

    // myfile2.close();
    // file2.close();
    // anotherfile2.close();

    // //=======================================================================================


    return 0;
}
