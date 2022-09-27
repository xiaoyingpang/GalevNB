#ifndef __PZ_CONFIG_H__
#define __PZ_CONFIG_H__

#include "filter.h"

using namespace std;


struct config_param
{
    std::string name;
    std::string value;
    std::string full_line;
};


#define NORMALIZATION_VONLY 1
#define NORMALIZATION_IONLY 2
#define NORMALIZATION_HONLY 3
#define NORMALIZATION_VH    4
#define NORMALIZATION_LBOL  20

class cconfig
{

 private:
    std::vector<config_param>  parameters;
    bool add_parameter(string name, string value, string full_line);
    
 public:
    virtual string get_parameter(string paramname);
    
    string isocont_file;
    string output_dir;

    string color_red;
    string color_blue;
    string magnitude;
    
    int nstar_max;
    int nstar_min;

    double scalingfactor;
    double timestep;
    double timestart;
    double timeend;
    
    double error_a;
    double error_b;
    double error_c;
    double error_d;
    
    double scale_col;
    double scale_mag;

    double imf_slope_low;
    double imf_slope_med;
    double imf_slope_hi;
    double imf_lowermass;
    double imf_uppermass;
    double imf_masslowmed;
    double imf_massmedhi;

    double noise_limit;

    string isolib;
    string speclib;
    double spectra_lambda_min;
    double spectra_lambda_max;

    cconfig();
    ~cconfig();
    int read_configfile(char* filename);
    int initialize();
    int checkconfig();

    string filtername_red;
    string filtername_blue;
    string filtername_mag;
    double zeropoint_red;
    double zeropoint_blue;
    double zeropoint_mag;

    filter* f_blue;
    filter* f_red;
    filter* f_mag;

    string filename_iso_spec;
    string filename_iso_int;

    int normalization_mode;
    bool integrate_by_weights;
    
};




#endif
