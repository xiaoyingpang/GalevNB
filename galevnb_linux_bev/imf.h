/* 
***********************************************************

  IMF class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/


#ifndef __IMF_H__
#define __IMF_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#define IMFMODE_STARS 1
#define IMFMODE_MASS  2

double integrate_powerlaw(double mass_low, double mass_hi, double slope);

class cimf
{


 protected:
    double lowermass;
    double uppermass;
    
    double slope_low;
    double slope_med;
    double slope_hi;

    double break_lowmed;
    double break_medhi;

    double scalingfactor_med;
    double scalingfactor_hi;
    double normalization;
    double slope_offset;

    virtual double compute_normalization();
    virtual double integrate_imf(double ml, double mu, double delta_slope);
    
 public:
    cimf();
    
    cimf(double lm, double um, double sl, double sm, double sh, double blm, double bmh);
    ~cimf();

    virtual double get_maxmass();
    virtual double get_minmass();
    virtual double integrate_imf_mass(double ml, double mu);
    virtual double integrate_imf_nstars(double ml, double mu);
};

// #include "imf_generalpowerlaw.h"

#endif
