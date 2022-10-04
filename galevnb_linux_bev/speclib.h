/* 
***********************************************************

  Spectra class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECLIB_H__
#define __SPECLIB_H__

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

#define SPECLIB_LEJEUNE 1
#define SPECLIB_UNKNOWN 0
#define SPECLIB_ASCIIPRECALIB 2
#define SPECLIB_BINARY 3

#define ILLEGAL_SPECID 0xffffffff

struct logg_bin
{
    double logg;
    spectrum* spec;
    double weight;
};

struct teff_bin
{
    double teff;
    vector< vector<logg_bin> > t;
};

struct metall_bin
{
    vector< vector<teff_bin> > m;
    double metallicity;
};

class spec_library
{
 protected:
    string index_filename;
    string lambda_filename;
    string lib_directory;
    int speclib_id;
    double lambda_min;
    double lambda_max;

    int skip_lambdas_front;
    int skip_lambdas_end;
    int read_lambdas;

    bool prepare_library();

    filter_library* filterlib;
    //vector<filter> filters;
    
    vector<spectrum*> library;
    vector<double> weights;
    vector<spectrum*> weights_additional_spectra;

    //vector< vector< vector<spectra*> > > slibrary;
    vector< vector<metall_bin> > slibrary;

    bool well_prepared;

    spectrum* get_spec(int spec_id);
    int initialize_library_struct();
    int get_bin_coords(double feh, double teff, double logg,
                       int* bin_m, int* bin_t, int* bin_g);

    int encode_library_pos(int bin_m, int bin_m2, int bin_t, int bin_t2,
                           int bin_g, int bin_g2);
    int decode_library_pos(int pos, int* bin_m, int* bin_m2, int* bin_t,
                           int* bin_t2, int* bin_g, int* bin_g2);

    double* weight(int index);
    
 public:
    vector<double> wavelength;
    spectrum* find_spec(
        double metall, double logg, double teff, double alpha_enh, 
        int s_metall,  int s_logg,  int s_teff,  int s_alpha_enh,
        int* index);


    spectrum* interpolate(spectrum* spec_low, spectrum* spec_hi,
                          double intfac, string normfilter);
    spectrum* interpolate(spectrum* spec_low, spectrum* spec_hi, string value,
                          double target, string normfilter, double* interfac);
    
    spec_library();
    spec_library(string galev_directory, string configfile, string speclib, double l_min,
                 double l_max, filter_library* filtlib);
    ~spec_library();
    bool is_valid();
    spectrum* get_spec(double metall, double logg, double teff,
                       double alpha_enh, string normfilter, double normmag);
    int add_to_library(spectrum* spec);
    spectrum* new_spec();

    bool reset_weights();
    spectrum* integrate_weights_to_spectrum();
    double get_spec_weights(double metall, double logg, double teff,
                            double alpha_enh, string normfilter,
                            double normmag, bool by_weight,
                            double additional_weight_factor);
    
    bool dump_structure();
    
    bool test();
    string info();
};

#endif

    

/* Structure of the slibrary vector ******************************
 *
 * metall_bin:                  < -2.0    index 0
 *               -2.0 <= [Fe/H] < -1.0          1
 *               -1.0 <=        <  0.0          2
 *                0.0 <=                        3
 *
 * ---------------------------------------------------------------
 *                
 * teff:         2000 <= Teff < 3000            0
 *               3000 <= Teff < 4000            1
 *               .
 *               .
 *               .
 *              49000 <= Teff < 50000           47
 *              50000 <= Teff ------------> Blackbody
 *
 * ----------------------------------------------------------------
 * 
 * logg                  logg < -2.0            0
 *               -2.0 <= logg < -1.0            1
 *               .
 *               .
 *                5.0 <= logg < 6.0             8
 *                6.0 <= logg                   9
 * 
 ***************************************************************** */
