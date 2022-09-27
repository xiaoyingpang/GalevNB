
#ifndef __LYCPHOTONS_H__
#define __LYCPHOTONS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "constants.h"

#define LUMCLASS_V   0
#define LUMCLASS_IV  1
#define LUMCLASS_III 2
#define LUMCLASS_II  3
#define LUMCLASS_I   4

using namespace std;

struct lyc_data
{
    double teff;
    double log_nlyc;
    double XXX;
    string spectral_type;
};

struct lyc_lumclass
{
    bool valid;
    vector<double> teff;
    vector<double> log_nlyc;
    gsl_interp_accel* acc;
    gsl_interp* g_interp;

    double teff_min;
    double teff_max;    
};


    
class lycphotons
{
  private:
    double photons;
    string datafile;
    bool loaded_and_valid;
    vector< vector<lyc_data> >data;
    vector< lyc_lumclass > lyc;

  public:
    lycphotons();
    ~lycphotons();
    bool load_data(string filename);
    
    double add_star(double Teff, double logg, double metallicity,
                    double log_bol, double weight);
    double get_integrated_flux();
    double reset_flux();
    
};



#endif
