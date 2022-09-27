/* 
***********************************************************

  Spec-Integrator

  written by: Ralf Kotulla
  (c) 2008-2010, The GALEV Team

  modifed by: Daniel Bialas 2011
  for the use in the Nbody-Programm

  modified by Christoph Euler 2012
  for use with the Giersz MC code

***********************************************************
*/

#include "cii.h"
#include <sys/time.h>
#include <sys/resource.h>
#include "speclib.h"

#ifndef SVN_REV
#define SVN_REV "???"
#endif


// Declare global objects.
int number_of_filters = 0;
const int array_size = 20;
int zeropoint_ids[array_size];
double zeropoint_mags[4][array_size];


struct star2mag_option
{
    vector<filter*> ugrizH;
    filter_library* filterlib;

    spec_library *speclib;
    
};

star2mag_option* opt = new star2mag_option;


star2mag_option* star2mag_initialize()
{
    star2mag_option* s2mopt = new star2mag_option;

    // Get path of source files via preprocessor macro defnition (in Makefile) and stringification (YEAH!).
#ifndef PATH_GALEV
  #error "PATH_GALEV hasn't been defined"
#endif
    // That's the stringification.
#define xstr(a) str(a)
#define str(a) #a
#define PATH_GALEV_STR xstr(PATH_GALEV)

    // Assign path to local variable.
    string path_galev = PATH_GALEV_STR;

    // Initialize object with standard filters.    
    s2mopt->filterlib = new filter_library();
    s2mopt->filterlib->add_filter(path_galev+"/standard_filters/Johnson_V.dat","V");
    s2mopt->filterlib->add_filter(path_galev+"/standard_filters/Cousins_Ic.dat","I");
    s2mopt->filterlib->add_filter(path_galev+"/standard_filters/BesselBrett_H.dat","H");

    // Initialize template spectra.    
    s2mopt->speclib = new spec_library(path_galev, "speclib.list", "lejeune_pc", 0, 1e6, s2mopt->filterlib);

    
    // Initialize user filter response curves.
    // Read configuration file 'filterlist.dat' and extract non-commented filter parameters.
    string filename = "filterlist.dat";
    
    // Accessing file.
    ifstream filterfile(filename.c_str(), ios_base::in);

    // Check file.
    if (!filterfile.good()) {
        cout << "Unable to open list of filters " << filename 
             << ", aborting ..." << endl;
    }

    string line;
    string token;

    int number_of_tokens;

    string filter_name[array_size];
    string data_file[array_size];
    
    while (getline(filterfile, line, '\n')) {
      if (line.size() <= 1) continue;
      if (line.at(0) == '#') continue;
      
      istringstream iss(line);
      iss >> filter_name[ number_of_filters ] >> data_file[ number_of_filters ] >> zeropoint_ids[ number_of_filters ] >> zeropoint_mags[ 0 ][ number_of_filters ] >> zeropoint_mags[ 1 ][ number_of_filters ] >> zeropoint_mags[ 2 ][ number_of_filters ] >> zeropoint_mags[ 3 ][ number_of_filters ];
      zeropoint_ids[ number_of_filters ]--; // Take into account that lowest ID is "1" while arrays begin with "0". 
      number_of_filters++;
    }
    filterfile.close();

    int entry = 0;

    s2mopt->ugrizH.resize(number_of_filters);
    while ( entry < number_of_filters ) {
      s2mopt->ugrizH.at( entry ) = new filter( path_galev + "/filter_response_curves/" + data_file[ entry ], filter_name[ entry ] );
      entry++;
    }

    return s2mopt;
}


// Pass number of filters.
void number_of_filters__pass( int *out__number_of_filters )
{
    // Pass globally stored number of filters.
    *out__number_of_filters = number_of_filters;
}


// calculate integrated magnitude from integrated spectrum of the cluster
void star2mag(star2mag_option* s2mopt,
              double* p_lteff, double* p_mstar, double* p_metallicity, double* p_lbol,
              double mags[array_size])
{
    double flux;

    double teff = *p_lteff,
        mstar = *p_mstar,
        metallicity = *p_metallicity,
    lbol = *p_lbol;
    
    spectrum* spec;

    double lteff=log10(teff);
    
    double loglbol = log10(lbol);

    // Compute logg from mstar and logr
    double logg = LOG_GRAVCONST
        + LOG_SUNMASS + log10(mstar)
        + LOG_SIGMA
        + (4 * lteff)
        + LOG_4PI 
        - (loglbol + LOG_SUNLUMBOL);

    // Convert Lbol into Mbol
    double mbol = -2.5 * loglbol + 4.72;// SUNMAG_BOL;

    //cout << "star2mag() further sends teff logg z" << teff << " " << logg << " " << metallicity << endl;

    s2mopt->speclib->reset_weights();
    s2mopt->speclib->get_spec_weights(metallicity, logg, teff, 0.0, "MBol", mbol, true, 1);
    spec = s2mopt->speclib->integrate_weights_to_spectrum();
    
    for (int i=0; i<number_of_filters; i++) {
        flux = spec->integrate_flux(s2mopt->ugrizH.at(i));
        mags[i] = -2.5 * log10(flux) - zeropoint_mags[ zeropoint_ids[ i ] ][ i ];
    }

    delete spec;
}
	

struct staropts
{
    double teff;
    double logg;
    double metallicity;
    double mbol;
    double weight;
};

struct specint_option{

    filter_library* filterlib;
    spec_library *speclib;

};


specint_option* opt_int = new specint_option; 


specint_option* specint_initialize(){
    specint_option* specintopt = new specint_option;

    // Get path of source files via preprocessor macro defnition (in Makefile) and stringification (YEAH!).
#ifndef PATH_GALEV
  #error "PATH_GALEV hasn't been defined"
#endif
    // That's the stringification.
#define xstr(a) str(a)
#define str(a) #a
#define PATH_GALEV_STR xstr(PATH_GALEV)

    // Assign path to local variable.
    string path_galev = PATH_GALEV_STR;

    // Set filter library for spectra.
    specintopt -> filterlib = new filter_library();
    specintopt -> filterlib -> add_filter(path_galev+"/standard_filters/Johnson_V.dat","V");
    specintopt -> filterlib -> add_filter(path_galev+"/standard_filters/Cousins_Ic.dat","I");
    specintopt -> filterlib -> add_filter(path_galev+"/standard_filters/BesselBrett_H.dat","H");

    specintopt->speclib = new spec_library(path_galev, "speclib.list", "lejeune_pc", 0, 1e6, specintopt->filterlib);

    return specintopt;
}

void reset_weights(specint_option* specintopt)
{
    specintopt->speclib->reset_weights();
}


/* add_star needs temp in nolog, mass in msun nolog, lum in nolog and z with 0.02 solar*/
void add_star(specint_option* specintopt,double lteff ,double mstar,double lbol,double Zmet)
{

    double teff, logg, mbol, metallicity;
    int    weight;
 
    teff = lteff;
    lteff = log10(lteff);
    double loglbol = log10(lbol);

	logg = LOG_GRAVCONST
               + LOG_SUNMASS + log10(mstar)
               + LOG_SIGMA
               + (4 * lteff)
               + LOG_4PI
               - (loglbol + LOG_SUNLUMBOL);

	  // Convert lbol into mbol
	 mbol = -2.5 * loglbol + 4.72; //SUNMAG_BOL;
  
	 weight = 1;

	 metallicity=Zmet;

	 specintopt->speclib->get_spec_weights(metallicity, logg, teff, 0.0,"MBol", mbol, true, weight);
}

void spec_output(specint_option* specintopt, string out)
{

    spectrum* spec;
    spec = specintopt->speclib->integrate_weights_to_spectrum();

    string output_spec(out);
    output_spec.erase(output_spec.find_last_not_of(" \n\r\t")+1); // Yohai's fix for trailing space in file names
    ofstream outf(output_spec.c_str(), ios_base::out);

    // Write header.
    outf << "#" << setw(11) << "wavelength" << setw(17) << "flux" << endl;
    outf << "#" << setw(11) << "[Angstrom]" << setw(17) << "[erg/s/cm^2/A]" << endl;

    //outf << "# lambda flux" << endl;
    for (int cl=0; cl<spec->extract_lambdas()->size(); cl++)
    {
        outf << scientific << setprecision(4) << setw(12) << spec->extract_lambdas()->at(cl)
             << scientific << setprecision(4) << setw(17) << spec->extract_fluxes()->at(cl) << endl;
    }
    outf.close();
    delete spec;

}

void spec2mag(specint_option* specintopt, star2mag_option* s2mopt, double mags[array_size] )
{
    double flux;

    spectrum* spec;
    spec = specintopt->speclib->integrate_weights_to_spectrum();
    
    for (int i=0; i<number_of_filters; i++) {
        flux = spec->integrate_flux(s2mopt->ugrizH.at(i));
        mags[i] = -2.5 * log10(flux) - zeropoint_mags[ zeropoint_ids[ i ] ][ i ];
    }

    delete spec;
}

extern "C" {

	void startomaginit_( int *out__number_of_filters )
        {
	     opt = star2mag_initialize();

	     // Assign number of filters used (for correct allocation in host routine).
	     number_of_filters__pass( out__number_of_filters );
	}

	void startomag_(double* p_lteff, double* p_mstar,double* p_lbol,double* p_Zmet,double* mag_array)
        {
             star2mag( opt, p_lteff, p_mstar, p_Zmet, p_lbol, mag_array);
	}

     void specint_initialize_()
     {
          opt_int = specint_initialize();
     }

     void reset_weights_()
     {
          reset_weights(opt_int);
     }

     void add_star_(double* p_lteff ,double* p_mstar,double* p_lbol,double* p_Zmet)
     {
          add_star(opt_int, *p_lteff , *p_mstar, *p_lbol, *p_Zmet);
     }

     void spec_output_(char *filename, int len_filename)
     {
          // Pass filename as parameter to spectra dump routine.
          spec_output(opt_int, filename);
     }


    void spec2mag_(double* mag_array)
     {
          spec2mag(opt_int,opt,mag_array);
     } 
}

