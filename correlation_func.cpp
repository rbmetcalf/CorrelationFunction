
/*
 * multiplane.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: bmetcalf & mpellejero
 *
 *      This program is for testing the adaptive griding on a random field by adapting to the
 *      high convergence regions.
 */

#include <slsimlib.h>
#include <standard.h>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include <time.h>

//#include <geometry.h>
#include "elliptic.h"
#include "gridmap.h"

using namespace std;

static std::mutex barrier;

//**************** test ************************

struct FunctorType1{
  FunctorType1(COSMOLOGY *cosmo,double radius,double redshift)
  : cosmology(cosmo),z(redshift),r(radius)
  {
    norm = 0.5/pi/pi;
  };
  
  COSMOLOGY *cosmology;
  double z;
  double r;
  
  double norm;
  
  double operator () (double lnk) {

    double k = exp(lnk);
    
    return norm*k*k*k*cosmology->powerCDMz(k,z);
  }
  
};

  


//***********************************************/


int main(int arg,char **argv){
  
  Utilities::print_date();
  
  COSMOLOGY cosmo;
  
  
}


