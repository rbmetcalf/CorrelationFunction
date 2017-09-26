
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

struct FunctorType{
  FunctorType(COSMOLOGY *cosmo,double radius,double redshift)
  : cosmology(cosmo),z(redshift),r(radius)
  {
    norm = 0.5/pi/pi/(1+z)/(1+z);
  };
  
  COSMOLOGY *cosmology;
  double z;
  double r;
  
  double norm;
  
  double operator () (double k) {
    
    double rk = r*k;
    double jo = sin(rk)/rk;
    if(rk < 1.0e-3) jo = 1.0;
    
    return norm*jo*k*k*cosmology->powerCDMz(k,z);
  }
};


struct CorrelationFunction{
  CorrelationFunction(COSMOLOGY *cosmo,double redshift
                      ,double k_max = 100,double k_min = 1.0e-3)
  : cosmology(cosmo),func(cosmo,1.0,redshift),k_max(k_max),k_min(k_min)
  {
    if(k_max < k_min) swap(k_min,k_max);
  };
  
  COSMOLOGY *cosmology;
  FunctorType func;
  double k_max;
  double k_min;
  
  double operator () (double radius) {
    double a = pi/radius/2;
    
    func.r = radius;
    
    double kmin = k_min;
    double kmax = MIN(a,k_max);
    double tmp,ans=0;
    do{
      tmp = Utilities::nintegrate<FunctorType,double>(func,kmin,kmax,1.0e-4);
      ans += tmp;
      kmin = kmax;
      kmax = MIN(kmax + a,k_max);
    }while(abs(tmp/ans) > 1.0e-4);
    
    return ans;
  }
};
  


//***********************************************/


int main(int arg,char **argv){
  
  Utilities::print_date();
  
  COSMOLOGY cosmo;
  double r = 1.3,z = 2;
  
  ofstream outputfile("correlation.csv");

  outputfile << "z,r,C" << std::endl;
  for(z=0.0;z <= 3;z += 1.0){
    CorrelationFunction correlation(&cosmo,z);

    for(r=0.001;r<5.0e2;r *=1.05){
      double corr = correlation(r);
    
      outputfile << z << "," << r << "," << corr << std::endl;
    
      std::cout << z << " , " << r << " , " << corr << std::endl;
    }
  }
  
  outputfile.close();
  
  std::cout << "finished" << std::endl;
}


