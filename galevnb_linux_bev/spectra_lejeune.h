/* 
***********************************************************

  Spectra class to be used as part of the GALEV package
  adapted for Lejeune spectra

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECTRA_LEJEUNE_H__
#define __SPECTRA_LEJEUNE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;

#include "filterlib.h"
#include "filter.h"
#include "spectra.h"


class spectrum_lejeune : public spectrum
{

 protected:
//    vector<double>* wavelengths;
//    vector<double> flux;
//    string filename;
//    bool   loaded;

//    double lambda_min;
//    double lambda_max;
//    int skip_lambdas_front;
//    int skip_lambdas_end;
//    int read_lambdas;
    bool read_specfile();

 public:
    spectrum_lejeune(string p_filename, vector<double> *wl,
             int *p_skip_begin, int *p_read, int *p_skip_end,
             filter_library* filtlib, void* sl);

    spectrum_lejeune(string p_filename, 
             double p_teff, double p_logg, double p_metall, double p_alpha_enh,
//           double p_lambda_min, double p_lambda_max,
             int p_skip_begin, int p_read, int p_skip_end, 
             vector<double>* p_wl, bool load_now,
             filter_library* filtlib, void* sl);

    void destroy();
    
//    double Teff;
//    double logg;
//    double metallicity;
//    double alpha_enh;

//    void dump();
};


#endif
