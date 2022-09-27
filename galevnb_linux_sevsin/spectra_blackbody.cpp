/* 
***********************************************************

  Spectra class to be used as part of the GALEV package
  adapted for Blackbody spectra

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#include "spectra_blackbody.h"
#include "constants.h"

bool spectrum_blackbody::read_specfile()
{
    double lambda_cm;
    
    if (!loaded) {
        // cout << "Computing BB with T=" << Teff << endl;

        unload(); // just to be sure we are not leaking memory
        flux = new vector<double>;
        flux->resize(wavelength->size());

        for (int i=0; i<flux->size(); i++) {
            lambda_cm = wavelength->at(i) * 1e-8;
            flux->at(i) = 4 * PI * PLANCKH * CVACUUM * CVACUUM * pow(lambda_cm,-5.0) / 
                //       -----------------------------------------------------------------
                         (exp((PLANCKH * CVACUUM)/(lambda_cm * KBOLTZMANN * Teff)) - 1.0);
        }
        
        loaded = true; 
        valid = loaded;

        for (int i=0; i<filterlib->get_number_filters(); i++) {
            set_normflux(filterlib->get_filter(i)->get_name(),
                         integrate_flux(filterlib->get_filter(i)));
        }
        set_normflux("MBol", integrate_flux());
    }
    
    return loaded;

}


spectrum_blackbody::spectrum_blackbody(double temp, vector<double> *wl,
                   filter_library *filtlib, void* sl)
{
    double new_lambda;
    string dummy, line;

    double lambda_min, lambda_max;

    wavelength = wl;
    filterlib = filtlib;
    speclib = sl;

    Teff = temp;
    logg = -99;
    metallicity = -99;
    alpha_enh = -99;

    flux = NULL;
    unload();
    
//    cout << "Creating blackbody spectrum with T=" << temp << endl;

    norm_fluxes.clear();
    norm_fluxes_names.clear();
    for (int i=0; i<norm_fluxes.size(); i++) {
        set_normflux(filterlib->get_filter(i)->get_name(), 0.0);
    }
    set_normflux("MBol", 0.0);

    loaded = false;   
    read_specfile();
}


