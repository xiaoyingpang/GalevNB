#include "spectra_blank.h"
#include "speclib.h"

spectrum_blank::spectrum_blank()
{
    //  cout << "Creating base class spectrum #2" << filename << "..." << endl;
    //  loaded = true;

    //
    // Nothing to do here, the work is done in the child class !
    //
//    spec_counter++;
//    cout << spec_counter << " ";
    flux = NULL;
    
    return;
}


spectrum_blank::~spectrum_blank()
{

    //    delete &flux;
    //    delete &norm_fluxes;
    //> 2013-08-13 C.O.
    // Fixing memory leak:
    // Unlike other spectral classes the class 'spectrum_blank' generates a new object 'wavelength' in addition to 'flux'.
    // The destructor 'unload' only takes care of the latter, so the following block of code is required to destroy the former.
    if (wavelength != NULL) {
        delete wavelength;
        wavelength = NULL;
    }
    //< 2013-08-13 C.O.

    unload();
    
    return;
}

spectrum_blank::spectrum_blank(spectrum* spec)
{
    // cout << "Creating blank spectrum ... " << flush;
    // if ((int)spec == 0) cout << "Found negative pointer " << endl;
    // cout << "Reading Teff = " << spec->Teff << endl;
    flux = NULL;
    unload();
        
    norm_fluxes.resize(spec->norm_fluxes.size());
    norm_fluxes_names.resize(spec->norm_fluxes_names.size());
    for (int i=0; i<norm_fluxes.size(); i++) {
        norm_fluxes.at(i) = spec->norm_fluxes.at(i);
        norm_fluxes_names.at(i) = spec->norm_fluxes_names.at(i);
    }

    wavelength = new vector<double>;
    wavelength->clear();
    flux = new vector<double>;
    
    Teff = spec->Teff;
    logg = spec->logg;
    metallicity = spec->metallicity;
    alpha_enh = spec->alpha_enh;
    valid = true;

    filename = "";
    loaded = true;

    // cout << "done !" << endl;
    
    return;
}


