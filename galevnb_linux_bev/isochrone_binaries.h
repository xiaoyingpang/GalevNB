/* 
***********************************************************

  Isochrone class to be used as part of the GALEV package
  inheritance to work with BINARIES-isochrones
  (cite Bertelli et al. 1994)

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __ISOCHRONE_BINARIES_H__
#define __ISOCHRONE_BINARIES_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "isochrone.h"

class isochrone_age_mass_binaries : public isochrone_age_mass 
{

 protected:
    double risoage;
    double smassm;
    double Teff;
    double logL;

    double weight;
    double logg;
    
    vector<double> magnitudes;

 public:
    isochrone_age_mass_binaries(string line);
    ~isochrone_age_mass_binaries();
    
    double get_mass_now();
    double get_mass_zams();
    double get_age();
    double get_logage();

    double get_Teff();
    double get_logg();
    double get_logL();

    double get_weight(cimf* imf, int weightmode,
                      double prevmass, double nextmass);
    
    double get_magnitude(string filtername);
};




#endif
