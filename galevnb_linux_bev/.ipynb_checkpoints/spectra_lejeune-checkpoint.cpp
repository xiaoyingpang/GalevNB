/* 
***********************************************************

  Spectra class to be used as part of the GALEV package
  adapted for Lejeune spectra

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#include "spectra_lejeune.h"


bool spectrum_lejeune::read_specfile()
{
  //  cout << "read_specfile " << filename << endl;

    if (!loaded) {
        cout << "Loading (Lejeune) " << filename << " ..." << endl << flush;
        std::ifstream specfile(filename.c_str(), std::ios_base::in);
        string line;

        if (!specfile.good()) {
            cerr << "Unable to open spec-file " << filename << endl;
            if (specfile.eof()) cout << "End of file reached" << endl;
            if (specfile.fail()) cout << "Operation failed" << endl;
            if (specfile.bad()) cout << "Operation is bad" << endl;
            return false;
        }
        // cout << "file opened ok" << endl;

        int curline = 0;
        
        for (int i=0; i<skip_lambdas_front; i++) {
            // Read as many lines as are not needed
            getline(specfile, line, '\n');
        }
        // cout << "skipped some lines" << endl;

        unload(); // just to be sure to clean memory
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

        // cout << "read a lot of flux values" << endl;

        for (int i=0; i<filterlib->get_number_filters(); i++) {
            // norm_fluxes.at(i) = integrate_flux(&filters->at(i));
            norm_fluxes.at(i) = integrate_flux(filterlib->get_filter(i));
            // cout << "Norm-flux " << i << " = " << norm_fluxes.at(i) << endl;
        }
        //        cout << "Integrated norm_fluxes" << endl;
        
        set_normflux("MBol", integrate_flux());
        //        cout << " done! (" << read_lambdas << ")" << endl;
    } else {
        //        cout << "File " << filename << " was already loaded ..." << endl;
    }

    //  cout << "read_specfile (Ende)" << endl;
    // cout << " done!" << endl;
    return loaded;

}

spectrum_lejeune::spectrum_lejeune(string p_filename, 
         double p_teff, double p_logg, double p_metall, double p_alpha_enh,
         int p_skip_front, int p_read, int p_skip_end, 
         vector<double>* p_wl, bool load_now,
         filter_library* filtlib, void* sl)
{

    // cout << "Creating lejeune spec " << flush;
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

    if (flux != NULL) {
        delete flux;
    }
    
//    filters = ff;
    filterlib = filtlib;

    loaded = false;

    // cout << "Init normfluxes ..." << endl;
    
    norm_fluxes.resize(filterlib->get_number_filters());
    norm_fluxes_names.resize(filterlib->get_number_filters());
    for (int i=0; i<filterlib->get_number_filters(); i++) {
        norm_fluxes.at(i) = 0;
        norm_fluxes_names.at(i) = filterlib->get_filter(i)->get_name();
        // cout << i << " " << filterlib->get_filter(i)->get_name();
    }
    norm_fluxes_names.push_back("MBol");
    norm_fluxes.push_back(0.0);
    if (load_now) loaded = read_specfile();

    // cout << "Lejeune done" << endl;
}


/*
void spectrum_lejeune::dump()
{
    cout << "Trying to dump (Lejeune) file ..." << endl;
//    read_specfile();
    if (read_specfile()) {
        cout << "Dumping " << flux.size() << " entries ..." << endl;
        for (int i=0; i<flux.size(); i++) {
            cout << i+1 << "  " 
                 << wavelength->at(i) << "  " 
                 << flux.at(i) << endl;
        }
        cout << "Dumping done!" << endl;
    }
    return;
}
*/

void spectrum_lejeune::destroy()
{

    return;
}


spectrum_lejeune::spectrum_lejeune(string p_filename, vector<double> *wl,
                   int *p_skip_begin, int *p_read, int *p_skip_end,
                   filter_library* filtlib, void* sl)
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


