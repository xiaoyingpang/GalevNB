/* 
***********************************************************

  Spectra class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECTRA_H__
#define __SPECTRA_H__

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

class spectrum
{
 protected:
    string filename;
    bool   loaded;
    vector<double>* wavelength;

//    double lambda_min;
//    double lambda_max;
    int skip_lambdas_front;
    int skip_lambdas_end;
    int read_lambdas;
    bool valid;
    virtual bool read_specfile();
    spectrum* add(spectrum* spec);
    void *speclib;
    
    filter_library* filterlib;
    vector<double>* flux;
    
    double bolometric_flux;

    int set_normflux(string filter, double value);

 public:
    vector<double> norm_fluxes;
    vector<string> norm_fluxes_names;
    spectrum();

    spectrum(string p_filename, vector<double> *wl,
             int *p_skip_begin, int *p_read, int *p_skip_end,
             filter_library* filtlib, void* sl);

    spectrum(string p_filename, 
             double p_teff, double p_logg, double p_metall, double p_alpha_enh,
//           double p_lambda_min, double p_lambda_max,
             int p_skip_begin, int p_read, int p_skip_end, 
             vector<double>* p_wl, bool load_now,
             filter_library* filtlib, void* sl);

    spectrum(vector<double>* p_wl, void* sl);


    //> 2013-08-13 C.O.
    // Fixing memory leak:
    // Replacing previous simple destructor '~spectrum();' by virtual destructor.
    // This is to make sure that the child destructor is called in addition to the parent's. 
    virtual ~spectrum();
    //< 2013-08-13 C.O.

    spectrum* operator+(spectrum* spec);
    spectrum* operator-(spectrum* spec);

    spectrum* operator*(double factor);
    spectrum* operator/(double factor);
    spectrum* mult(double factor);
    void operator=(spectrum spec);

    double operator*(filter ff);
    double integrate_flux();
    double integrate_flux(filter *ff);
    double integrate_flux(double start, double end);
    void operator*=(const double factor);
    virtual void dump(int steps);
    virtual string info();
    
    bool is_valid();
    double get_normflux(string filter);
    double get_normflux(int filterid);

    virtual void destroy();
    virtual bool unload();
        
    double Teff;
    double logg;
    double metallicity;
    double alpha_enh;

    vector<double>* extract_lambdas();
    vector<double>* extract_fluxes();

};

#include "spectra_lejeune.h"
#include "spectra_blackbody.h"
#include "spectra_asciiprecalib.h"
#include "spectra_blank.h"
#include "spectra_binary.h"

#endif
