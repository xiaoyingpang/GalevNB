/* 
***********************************************************

  Spectra class to be used as part of the GALEV package
  adapted for Lejeune-like spectra, but with the integrated
  fluxes for the calibration filters read from index file

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#include "spectra_asciiprecalib.h"


bool spectrum_asciiprecalib::read_specfile()
{
    if (!loaded) {
        std::ifstream specfile(filename.c_str(), std::ios_base::in);
        string line;

        if (!specfile.good()) {
            cerr << "Unable to open spec-file " << filename << endl;
            if (specfile.eof()) cout << "End of file reached" << endl;
            if (specfile.fail()) cout << "Operation failed" << endl;
            if (specfile.bad()) cout << "Operation is bad" << endl;
            return false;
        }

        int curline = 0;
        
        for (int i=0; i<skip_lambdas_front; i++) {
            // Read as many lines as are not needed
            getline(specfile, line, '\n');
        }

        unload();
        flux = new vector<double>;
        flux->resize(read_lambdas);
        while (getline(specfile, line, '\n') && curline < read_lambdas) {
            if (line.size() <= 0) continue;
            if (line.at(0) == '#') continue;
            
            istringstream iss(line);
            iss >> flux->at(curline++);
        }
        specfile.close();

        loaded = (curline == read_lambdas);
        valid = loaded;
    }

    return loaded;
}

spectrum_asciiprecalib::spectrum_asciiprecalib(string p_filename, 
         double p_teff, double p_logg, double p_metall, double p_alpha_enh,
         int p_skip_front, int p_read, int p_skip_end, 
         vector<double>* p_wl, bool load_now,
         filter_library* filtlib, void* sl, string additional_parameters)
{

    // cout << "Creating asciiprecalib spectra from file " << filename << endl << flush;
    
    //    cout << "Creating lejeune spec " << flush;
    filename = p_filename;
    Teff = p_teff;
    logg = p_logg;
    metallicity = p_metall;
    alpha_enh = p_alpha_enh;
    wavelength = p_wl;
    speclib = sl;

    skip_lambdas_front = p_skip_front;
    skip_lambdas_end   = p_skip_end;
    read_lambdas       = p_read;

    flux = NULL;
    unload();
    flux = new vector<double>;

    filterlib = filtlib;

    loaded = false;

    // Deal with pre-computed norm-fluxes
    string dummy, filtername, entry;
    double log_flux;

    // First, replace all = by spaces in order to feed it to istringstream
    dummy = additional_parameters;
    for (int i=0;i<dummy.length(); i++)
        dummy.at(i) = dummy.at(i) == '=' ? ' ' : dummy.at(i);
    istringstream iss(dummy);

    // Now reset all values that might have been there (actually should be
    // empty already, just to make sure)
    norm_fluxes.clear();
    norm_fluxes_names.clear();
    // and finally, get filternames, fluxes and store into array
    while (!iss.eof()) {
        iss >> filtername >> log_flux;
        set_normflux(filtername, pow(10.0, -0.4*log_flux));
    }
    
    if (load_now) loaded = read_specfile();
}

void spectrum_asciiprecalib::destroy()
{

    return;
}


spectrum_asciiprecalib::spectrum_asciiprecalib(string p_filename, vector<double> *wl,
                   int *p_skip_begin, int *p_read, int *p_skip_end,
                   filter_library* filtlib, void* sl, string additional_parameters)
{
    double new_lambda;
    string dummy, line;

    double lambda_min, lambda_max;

    filterlib = filtlib;
    speclib = sl;

    skip_lambdas_front = skip_lambdas_end = read_lambdas = 0;
    std::ifstream lambdafile(p_filename.c_str(), std::ios_base::in);
    if (!lambdafile.good()) {
        cerr << "Unable to load speclib lambdafile " << p_filename << endl;
	//        exit(-1);
    }
    while(getline(lambdafile, line, '\n')) {
        if (line.size() <= 0) continue;
        if (line.at(0) == '#') continue;

        istringstream iss(line);
        iss >> dummy >> new_lambda;
        if (new_lambda < lambda_min) {
            skip_lambdas_front++;
        } else if (new_lambda >= lambda_min &&
                   new_lambda <= lambda_max) {
            read_lambdas++;
            wl->push_back(new_lambda);
        } else {
            skip_lambdas_end++;
        }
    }
    lambdafile.close();

    *p_skip_begin = skip_lambdas_front;
    *p_read = read_lambdas;
    *p_skip_end = skip_lambdas_end;
}


