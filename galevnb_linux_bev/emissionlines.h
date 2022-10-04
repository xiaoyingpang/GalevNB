

#ifndef __EMISSIONLINES_H__
#define __EMISSIONLINES_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>
#include <gsl/gsl_sf_erf.h>

#include "constants.h"
#include "spectra.h"
#include "lycphotons.h"

using namespace std;

struct linestrength
{
    vector<double> metallicity;
    vector<double> rel_strength;
    double wavelength;
    string name;
};

class emissionlines
{
  private:
    string definition_file;
    vector<linestrength> lines;

    double linewidth_fwhm;
    double flux_to_intensity;
    double line_sigma;
    
  public:
    emissionlines();
    ~emissionlines();

    int load_definitions(string filename);
    int add_to_spec(spectrum* spec, lycphotons* lyc, double metallicity);
};



#endif
