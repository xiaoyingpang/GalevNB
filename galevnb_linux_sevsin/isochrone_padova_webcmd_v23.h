/* 
***********************************************************

  Isochrone class to be used as part of the GALEV package
  inheritance to work with PADOVA-isochrones
  (cite Bertelli et al. 1994)

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __ISOCHRONE_PADOVA_WEBCMD_V23_H__
#define __ISOCHRONE_PADOVA_WEBCMD_V23_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "isochrone.h"

class isochrone_age_mass_padova_webcmd_v23 : public isochrone_age_mass
{

 protected:
    double log_age;
    double mass_initial;
    double mass_actual;
    double logL;
    double logT;
    double logG;
    double mbol;
    double c_o;
    double M_hec;

    double period;
    double pmode;
    double logMdot;
    double int_IMF;
    double stage;
    
    vector<double> magnitudes;

 public:
    isochrone_age_mass_padova_webcmd_v23();
    isochrone_age_mass_padova_webcmd_v23(string line);
    ~isochrone_age_mass_padova_webcmd_v23();
    
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
