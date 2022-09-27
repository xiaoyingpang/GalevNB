/* 
***********************************************************

  Isochrone class to be used as part of the GALEV package
  inheritance to work with PADOVA-isochrones
  (cite Bertelli et al. 1994)

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __ISOCHRONE_PADOVA_H__
#define __ISOCHRONE_PADOVA_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "isochrone.h"

class isochrone_age_mass_padova : public isochrone_age_mass 
{

 protected:
    double risoage;
    double smassm;
    double logT;
    double logL;

    double rflum;
    double rmwind;
    vector<double> magnitudes;

 public:
//    isochrone_age_mass_padova(string line);
    isochrone_age_mass_padova(string line);
    ~isochrone_age_mass_padova();
    
    double get_mass_now();
    double get_mass_zams();
    double get_age();
    double get_logage();

    double get_Teff();
    double get_logg();
    double get_logL();

    double get_magnitude(string filtername);
};




#endif
