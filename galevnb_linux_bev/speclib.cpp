#include "speclib.h"

string trim_string(string input)
{
    size_t first = input.find_first_not_of(" ");
    size_t last = input.find_last_not_of(" ");
    first = (first == string::npos) ? 0 : first;
    last = (last == string::npos) ? input.length() : last;
    
    return input.substr(first, last-first+1);
}

spec_library::spec_library()
{

}



spec_library::spec_library(string galev_directory, string listfile, string speclib, double l_min, double l_max, filter_library* filtlib)
{
    listfile =  galev_directory + "/" + listfile;

    std::ifstream libfile(listfile.c_str(), std::ios_base::in);
    string line;
    bool loaded_ok = false;

    if (!libfile.good()) {
        cerr << "Unable to open speclist-file " << listfile << endl;
    }

    lambda_min = l_min;
    lambda_max = l_max;

    string isofilename, dummy, readspeclib;
    double metallicity, alpha_enh;
    bool found_section = false;

    filterlib = filtlib;
    weights_additional_spectra.clear();
    
    //cout << "Reading speclib-overview ..." << flush;
    while (getline(libfile, line, '\n')) {
        if (line.size() <= 0) continue;
        if (line.at(0) == '#') continue;

        if (line.at(0) == '@') {
            istringstream iss(line);
            iss >> dummy >> readspeclib >> speclib_id;
            
            if (readspeclib == speclib) {
                found_section = true;
                { // Read directory
                    getline(libfile, line, '\n');
                    istringstream iss(line);
                    iss >> dummy;
		    lib_directory = galev_directory + "/" + dummy;
                }{ // Read filename with index
                    getline(libfile, line, '\n');
                    istringstream iss(line);
                    iss >> dummy;
                    index_filename = lib_directory + "/" + dummy;
                }{ // Read name of file containing wavelength points
                    getline(libfile, line, '\n');
                    istringstream iss(line);
                    iss >> dummy;
                    lambda_filename = lib_directory + "/" + dummy;
                }
                
                // cout << endl << "Lib-Dir. = " << lib_directory << endl;
                loaded_ok = prepare_library();
                break;
            }
        }
    }

    if (!found_section) {
        cout << "Unable to find specified library" << endl;
    }

    well_prepared = (found_section && loaded_ok);

//    cout << "Lib-Dir. = " << lib_directory << endl;
}

bool spec_library::prepare_library()
{

    // cout << "Preparing library ... " << endl;

    double Teff;
    double logg;
    double metallicity;
    double alpha_enh;
    string filename, dummy, additional_parameters;
    
    spectrum *spec_base;
    spectrum *spec_lejeune;
//    spectrum_lejeune *spec_lejeune;

    string line;

    /*
    filter johnsonV("input_data/filter/johnsonV.dat","V");
    filter cousinsI("input_data/filter/Cousins_Ic.dat","I");
    filter bessellH("input_data/filter/bessellH.dat","H");
    filters.clear();
    filters.push_back(johnsonV);
    filters.push_back(cousinsI);
    filters.push_back(bessellH);
    */
    initialize_library_struct();
    
    double new_lambda;
    skip_lambdas_front = skip_lambdas_end = read_lambdas = 0;
    std::ifstream lambdafile(lambda_filename.c_str(), std::ios_base::in);
    if (!lambdafile.good()) {
        cerr << "Unable to load speclib lambdafile " << lambda_filename << endl;
	//        exit(-1);
        return false;
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
            wavelength.push_back(new_lambda);
        } else {
            skip_lambdas_end++;
        }
        //wavelength.push_back(new_lambda);

    }
    lambdafile.close();
/*    for (int i=0; i<wavelength.size(); i++) {
        cout << wavelength.at(i) << "  ";
    }
    cout << endl << endl << endl;*/

    // cout << "Lib-Dir. = " << lib_directory << endl;
    std::ifstream index(index_filename.c_str(), std::ios_base::in);
    if (!index.good()) {
        cerr << "Unable to load speclib indexfile " << index_filename << endl;
	//        exit(-1);
    }
    
    while(getline(index, line, '\n')) {
        if (line.size() <= 0) continue;
        if (line.at(0) == '#') continue;

        istringstream iss(line);
        iss >> Teff >> logg >> metallicity >> alpha_enh >> dummy;
        filename = lib_directory + "/" + dummy;

        additional_parameters = iss.str().substr(iss.tellg());
        
        // cout << "Adding spectrum " << filename << " to library ..." << endl;
        switch (speclib_id) {
            case SPECLIB_LEJEUNE:
              // cout << "Reading Lejeune spectra" << endl;
                spec_lejeune = new 
                    spectrum_lejeune(filename, Teff, logg, metallicity, alpha_enh,
                                     skip_lambdas_front, read_lambdas, skip_lambdas_end, 
                                     &wavelength, false, filterlib, this);
                //spec_lejeune->dump();
                add_to_library(spec_lejeune);
                // cout << "Lejeune done" << endl;
                break;

            case SPECLIB_ASCIIPRECALIB:
                spec_lejeune = new 
                    spectrum_asciiprecalib(filename, Teff, logg, metallicity,
                                           alpha_enh, skip_lambdas_front,
                                           read_lambdas, skip_lambdas_end, 
                                           &wavelength, false, filterlib, this,
                                           trim_string(additional_parameters));
                add_to_library(spec_lejeune);
                break;
                
            case SPECLIB_UNKNOWN:
            default:
              // cout << "I'm here" << flush;
                spec_base = new 
                    spectrum(filename, Teff, logg, metallicity, alpha_enh,
                             skip_lambdas_front, read_lambdas, skip_lambdas_end, 
                             &wavelength, false, filterlib, this);
                add_to_library(spec_base);
                break;
        }
    }
    index.close();
    return true;   
}

int spec_library::get_bin_coords(double feh, double teff, double logg,
                                 int* bin_m, int* bin_t, int* bin_g)
{
    *bin_m = (feh < -2.0) ? 0 : ( (feh >= 0) ? 3 : floor(feh + 2.)+1);
    *bin_g = (logg < -2.0) ? 0 : (logg >= 6.0 ? 9 : floor(logg + 2)+1);
    *bin_t = (teff < 2000) ? 0 : (teff >= 50000 ? 48 : floor((teff-2000)/1000)+1);
    
    return 0;
}

int spec_library::add_to_library(spectrum* spec)
{
    double feh, logg, teff;
    int bin_m, bin_t, bin_g, bin_m2, bin_t2, bin_g2;
    bool found, insert;

    vector< vector<teff_bin> > *select_m;
    vector< vector<logg_bin> > *select_t;
    
    // Get parameters from spectrum
    feh = spec->metallicity;
    logg = spec->logg;
    teff = spec->Teff;

    // Now find out which pre-sort bins we need
    get_bin_coords(feh, teff, logg, &bin_m, &bin_t, &bin_g);
    // cout << feh << " " << teff << " " << logg << " --> "
    //      << bin_m << " " << bin_t << " " << bin_g << endl;
    
    // Check if this metallicity already has its branch
    found = insert = false;
    for (bin_m2 = 0; bin_m2<slibrary.at(bin_m).size(); bin_m2++) {
        if (slibrary.at(bin_m).at(bin_m2).metallicity == feh) {
            found = true;
            select_m = &slibrary.at(bin_m).at(bin_m2).m;
            break;
        }
        if (slibrary.at(bin_m).at(bin_m2).metallicity > feh) {
            insert = true;
            break;
        }
    }
    
    // Ok, we haven't found any entry with this metallicity yet, so let's create
    // a new branch, including the necessary Teff-bin subbranches
    if (!found) {
        // cout << "Creating new metall...." << flush;
        if (!insert) {
            // No larger metallicity available yet, so simple add
            // to the end of the list
            slibrary.at(bin_m).resize(slibrary.at(bin_m).size()+1);
            bin_m2 = slibrary.at(bin_m).size() - 1;
        } else {
            // There is already a higher metallicity, so insert into the
            // right position instead of adding to the end
            vector<metall_bin>::iterator it;
            it = slibrary.at(bin_m).begin();
            metall_bin mb;
            slibrary.at(bin_m).insert(it+bin_m2, mb);
        }
        select_m = &slibrary.at(bin_m).at(bin_m2).m;
        slibrary.at(bin_m).at(bin_m2).metallicity = feh;
        select_m->resize(49);
        for (int i=0; i<select_m->size(); i++) {
            select_m->at(i).clear();
        }
        // cout << "done !" << endl;
    }

    // Now repeat above game to find the right teff-bin, and create the
    // appropriate logg sub-branches if necessary
    found = insert = false;
    for (bin_t2 = 0; bin_t2 < select_m->at(bin_t).size(); bin_t2++) {
        if (select_m->at(bin_t).at(bin_t2).teff == teff) {
            found = true;
            select_t = &select_m->at(bin_t).at(bin_t2).t;
            break;
        }
        if (select_m->at(bin_t).at(bin_t2).teff > teff) {
            insert = true;
            break;
        }
    }
    if (!found) {
        // cout << "Creating new Teff=" << teff << "   M=" << bin_m << "."
        //      << bin_m2 << "   T=" << bin_t << endl;
        // As above with the metallicities, either push_back to create new
        // entry at the end of the list or insert into the middle
        if (!insert) {
            select_m->at(bin_t).resize(select_m->at(bin_t).size()+1);
            bin_t2 = select_m->at(bin_t).size() - 1;
        } else {
            // There is already a higher metallicity, so insert into the
            // right position instead of adding to the end
            vector<teff_bin>::iterator it;
            it = select_m->at(bin_t).begin();
            teff_bin tb;
            select_m->at(bin_t).insert(it+bin_t2, tb);
        }
        // Enter the value of the new temperature ...
        select_m->at(bin_t).at(bin_t2).teff = teff;
        // ... and create the sub-vectors for the logg's
        select_t = &select_m->at(bin_t).at(bin_t2).t;
        select_t->resize(10);
        for (int i=0; i<select_t->size(); i++) {
            select_t->at(i).clear();
        }
    }

    // Almost final step: Look for matching logg
    // Still to be decided: What to do in case the exact logg-value
    // already exists? Overwrite, skip, simply ignore, or add duplicate?
    found = insert = false;
    for (bin_g2 = 0; bin_g2 < select_t->at(bin_g).size(); bin_g2++) {
        if (select_t->at(bin_g).at(bin_g2).logg == logg) {
            found = true;
            // Now we found two spectra with the same parameters
            // (metallicity, Teff, logg), so skip this one.
            return -1; // <-- Change the return value
        }
        if (select_t->at(bin_g).at(bin_g2).logg > logg) {
            insert = true;
            break;
        }
    }
    if (!found) {
        // cout << "Creating new logg=" << logg << endl;
        // As above with the metallicities and Teff, either push_back to create new
        // entry at the end of the list or insert into the middle
        if (!insert) {
            select_t->at(bin_g).resize(select_t->at(bin_g).size()+1);
            bin_g2 = select_t->at(bin_g).size() - 1;
        } else {
            // There is already a higher logg, so insert into the
            // right position instead of adding to the end
            vector<logg_bin>::iterator it;
            it = select_t->at(bin_g).begin();
            logg_bin gb;
            select_t->at(bin_g).insert(it+bin_g2, gb);
        }

        // Now that all vector entries are created, enter the
        // logg value and the pointer to the spectrum
        select_t->at(bin_g).at(bin_g2).logg = logg;
        select_t->at(bin_g).at(bin_g2).spec = spec;
        select_t->at(bin_g).at(bin_g2).weight = 0;
    }

    // That's it, the spectrum is saved at its optimum position, let's
    // return the pointer to it that encodes its position in the library.

    return encode_library_pos(bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2);
}

spec_library::~spec_library()
{

}
spectrum* spec_library::new_spec()
{
    spectrum* newspec = new spectrum(&wavelength, this);
    return newspec;
}

inline int sign(double x) 
{
    if (x == 0) return 1;
    return (int)(x / fabs(x));
}


spectrum* spec_library::find_spec(double metall, double logg, double teff, double alpha_enh, 
                                  int s_metall,  int s_logg,  int s_teff,  int s_alpha_enh,
                                  int* index)
{
    // First find out which pre-sort bins we need
    int bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2;
    get_bin_coords(metall, teff, logg, &bin_m, &bin_t, &bin_g);

    bool vm_metall, vm_teff, vm_logg;  // vm = valid match
    int startbin, curbin;

    vector< vector<teff_bin> > *select_m;
    vector< vector<logg_bin> > *select_t;

    // For now set index to illegal value. If we find a match we will
    // set the correct value at the end.
    if (index != NULL) *index = ILLEGAL_SPECID;
    
    vm_metall = vm_teff = vm_logg = false;

    // cout << "Looking for star - [Fe/H]=" << metall << ", Teff="<< teff << ", logg=" << logg << endl;
    
    // Find the best-matching metallicity-index
    // cout << "Searching for metallicities ..." << flush;
    while (!vm_metall && bin_m >= 0 && bin_m < slibrary.size()) {
        vm_metall = false;
        startbin = (s_metall > 0) ? 0 : slibrary.at(bin_m).size()-1;
        for (bin_m2 = startbin; bin_m2 < slibrary.at(bin_m).size() && bin_m2>=0 ; bin_m2+=s_metall) {
            if ((s_metall * slibrary.at(bin_m).at(bin_m2).metallicity) >= (metall * s_metall)) {
                vm_metall = true;
                break;
            }
        }
        // If no matching spectrum found yet, then proceed to next bin
        if (!vm_metall) {
            bin_m += s_metall;
        }
    }
    if (!vm_metall) {
        // cout << "Not found" << endl;
        return NULL;
    }
    // cout << " done!" << endl;
    
    select_m = &slibrary.at(bin_m).at(bin_m2).m;
    // cout << "Searching for Teff ..." << flush;
    // Find the best-matching teff-index
    while (!vm_teff && bin_t >= 0 && bin_t < select_m->size()) {
        vm_teff = false;
        startbin = (s_teff > 0) ? 0 : select_m->at(bin_t).size()-1;
        for (bin_t2 = startbin; bin_t2 < select_m->at(bin_t).size() && bin_t2 >= 0; bin_t2+=s_teff) {
            if ((s_teff * select_m->at(bin_t).at(bin_t2).teff) >= (teff * s_teff)) {
                vm_teff = true;
                break;
            }
        }
        // If no matching spectrum found yet, then proceed to next bin
        if (!vm_teff) {
            bin_t += s_teff;
        }
    }
    if (!vm_teff) {
        // cout << "Not found" << endl;
        return NULL;
    }
    // cout << " done!" << endl;

    select_t = &select_m->at(bin_t).at(bin_t2).t;
    // cout << "Searching for logg ..." << flush;
    // Find the best-matching logg-index
    while (!vm_logg && bin_g >= 0 && bin_g < select_t->size()) {
        vm_logg = false;
        startbin = (s_logg > 0) ? 0 : select_t->at(bin_g).size()-1;
        for (bin_g2 = startbin; bin_g2 < select_t->at(bin_g).size() && bin_g2 >= 0; bin_g2+=s_logg) {
            if ((s_logg * select_t->at(bin_g).at(bin_g2).logg) >= (logg * s_logg)) {
                vm_logg = true;
                break;
            }
        }
        // If no matching spectrum found yet, then proceed to next bin
        if (!vm_logg) {
            bin_g += s_logg;
        }
    }
    if (!vm_logg) {
        // cout << "Not found" << endl;
        return NULL;
    }
    // cout << " done!" << endl;

    if (index != NULL) {
        *index = encode_library_pos(bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2);
    }
    
    return select_t->at(bin_g).at(bin_g2).spec;
    
/*
    if (index != NULL) *index = best_fit_index;
    
    if (one_found) {
        //cout << "Nearest match: Teff=" << best_fit->Teff << ", logg=" << best_fit->logg 
        //     << ", metall=" << best_fit->metallicity << " !" << endl << endl;
        return best_fit;
    }

    // cout << "Nothing found" << endl << endl;
    return best_fit;
*/
}

spectrum* spec_library::interpolate(spectrum* spec_low, spectrum* spec_hi, string value,
                                    double target, string normfilter, double* interfac)
{
    double diff, interpol_fac;
    
    if (spec_low == 0 && spec_hi == 0) {
        if (interfac != NULL) *interfac = -1;
        return 0; 
    }
    if (spec_hi == 0) {
        if (interfac != NULL) *interfac = 0;
        return *spec_low * 1.0;
    }
    if (spec_low == 0) {
        if (interfac != NULL) *interfac = 1;
        return *spec_hi * 1.0;
    }

    if (value == "Teff") {
        diff = log10(spec_hi->Teff) - log10(spec_low->Teff);
        interpol_fac = (log10(target) - log10(spec_low->Teff)) / diff;
    } else if (value == "logg") {
        diff = spec_hi->logg  - spec_low->logg;
        interpol_fac = (target - spec_low->logg) / diff;
    } else if (value == "metallicity") {
        // cout << "xxx " << spec_hi->metallicity << " " << spec_low->metallicity << endl;
        diff = spec_hi->metallicity  - spec_low->metallicity;
        // cout << "xxx " << target << " " << diff << endl;
        interpol_fac = (target - spec_low->metallicity) / diff;
        // cout << "xxx " << interpol_fac << endl;
    } else {
        cout << "Illegal option for interpolation";
        return 0;
    }
    
        
    if (diff == 0) {
        if (interfac != NULL) *interfac = 1;
        return *spec_hi * 1;
    } else {
        if (interpol_fac > 1) {
            if (interfac != NULL) *interfac = 1;
            return *spec_hi * 1;
        } else if (interpol_fac < 0) {
            if (interfac != NULL) *interfac = 0;
            return *spec_low * 1;
        } else {
            if (interfac != NULL) *interfac = interpol_fac;
            return interpolate(spec_low, spec_hi, interpol_fac, normfilter);
        }
    }
    
    return NULL;
}

spectrum* spec_library::interpolate(spectrum* spec_low, spectrum* spec_hi, double intfac, string normfilter)
{
    spectrum* result;
    if (spec_low == 0 && spec_hi == 0) return 0;
    if (spec_low == 0) spec_low = spec_hi;
    if (spec_hi  == 0) spec_hi = spec_low;

    double normfac_hi, normfac_low, normfac_result;
    double logfac_hi, logfac_low, logfac_result;

    double target_logg, target_Teff;

    // cout << "Computing interpol-facs ..." << endl;
    
    target_logg = (intfac  * spec_hi->logg) + ((1 - intfac) * spec_low->logg);
    target_Teff = (intfac  * spec_hi->Teff) + ((1 - intfac) * spec_low->Teff);

    logfac_low    = 4 * log10(spec_low->Teff) - spec_low->logg;
    logfac_hi     = 4 * log10(spec_hi->Teff)  - spec_hi->logg;
    logfac_result = 4 * log10(target_Teff)    - target_logg;

    normfac_low   = 1; //spec_low->get_normflux(normfilter); //pow(10.0, (logfac_result - logfac_low));
    normfac_hi    = 1; //spec_hi->get_normflux(normfilter); //pow(10.0, (logfac_result - logfac_hi));

    // spectrum* part1 = *spec_hi  * (     intfac / spec_hi->get_normflux(normfilter));
    // spectrum* part2 = *spec_low * ((1 - intfac)/spec_low->get_normflux(normfilter));

    // cout << "Computing new specs " << endl;
    if (intfac > 1 || intfac < 0)
        cout << "Interpolating with negative factor " << intfac << endl;
    
    spectrum* part1 = *spec_hi  * (     intfac / normfac_hi);
    spectrum* part2 = *spec_low * ((1.0 - intfac)/ normfac_low);

    result = *part1 + part2;
    part1->destroy();
    part2->destroy();
    
    result->Teff = target_Teff;
    result->metallicity 
        =      intfac  * spec_hi->metallicity
        + (1 - intfac) * spec_low->metallicity;

    result->logg = target_logg;

    return result;
}

spectrum* spec_library::get_spec(double metall, double logg, double teff, double alpha_enh,
                                 string normfilter, double normmag)
{
    int normflux_id = 0;
    spectrum* retspec;
    double target_flux = normmag;

    target_flux = pow(10.0, -0.4*normmag);

    if (teff > 50000) {
        spectrum* bbspec =
            new spectrum_blackbody(teff,&wavelength,filterlib,this);
        retspec = (*bbspec * (target_flux / bbspec->get_normflux(normfilter)));
        delete bbspec;
        return retspec;
    }

    double interpol_fac = 0, diff;
    spectrum lowT_lowM, lowT_hiM, hiT_lowM, hiT_hiM;
    spectrum lowM, hiM;
    spectrum *result;

    //
    // Find surrounding 8 best spectra 
    // Closest in logG, Teff and metall
    //

    vector<spectrum*> level1(8);
    vector<spectrum*> level2(4);
    vector<spectrum*> level3(2);
    vector<spectrum*> level4(1);
    
    level1.at(0) = find_spec(metall, logg, teff, alpha_enh, -1, -1, -1, +1, NULL);
    level1.at(1) = find_spec(metall, logg, teff, alpha_enh, -1, -1, +1, +1, NULL);
    level1.at(2) = find_spec(metall, logg, teff, alpha_enh, -1, +1, -1, +1, NULL);
    level1.at(3) = find_spec(metall, logg, teff, alpha_enh, -1, +1, +1, +1, NULL);
    level1.at(4) = find_spec(metall, logg, teff, alpha_enh, +1, -1, -1, +1, NULL);
    level1.at(5) = find_spec(metall, logg, teff, alpha_enh, +1, -1, +1, +1, NULL);
    level1.at(6) = find_spec(metall, logg, teff, alpha_enh, +1, +1, -1, +1, NULL);
    level1.at(7) = find_spec(metall, logg, teff, alpha_enh, +1, +1, +1, +1, NULL);

    //
    // for all 4 blocks
    // first interpolate in Teff
    // 
    for (int i=0; i<4; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level2.at(i) = interpolate(level1.at(low), level1.at(hi), "Teff", teff, normfilter, NULL);
    }

    //
    // for all 2 blocks
    // next interpolate in logg
    // 
    for (int i=0; i<2; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level3.at(i) = interpolate(level2.at(low), level2.at(hi), "logg", logg, normfilter, NULL);
    }
    for (int i=0; i<4; i++) {
        if (level2.at(i) != 0) delete level2.at(i);
    }
    
    //
    // for all 2 blocks
    // next interpolate in metallicity
    // 
    for (int i=0; i<1; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level4.at(i) = interpolate(level3.at(low), level3.at(hi), "metallicity", metall, normfilter, NULL);
    }
    for (int i=0; i<2; i++) {
        if (level3.at(i) != 0) delete level3.at(i);
    }

    retspec = (*level4.at(0) * (target_flux / level4.at(0)->get_normflux(normfilter)));
    if (level4.at(0) != 0) delete level4.at(0);
    
    return retspec;
}


double spec_library::get_spec_weights(double metall, double logg, double teff, double alpha_enh,
                                      string normfilter, double normmag, bool by_weight,
                                      double additional_weight_factor)
{
    int normflux_id = 0;
    spectrum* retspec;
    double target_flux = normmag;

    target_flux = pow(10.0, -0.4*normmag);

    // if (logg < -0.7) return 0;
    
    // cout << scientific << setprecision(5) << metall << " " << logg << " " << teff << " "
    //      << target_flux << " " << additional_weight_factor << endl;
    
    if (teff > 50000) {
        // cout << "Creating new BB spectrum for Teff=" << teff << flush;

        spectrum* bbspec =
            new spectrum_blackbody(teff,&wavelength,filterlib,this);

        // cout << "target=" << target_flux << "   " << normfilter
        //      << "=" <<  bbspec->get_normflux(normfilter) << endl;
        
        retspec = (*bbspec * (additional_weight_factor * 1e-10 *
                              target_flux / bbspec->get_normflux(normfilter)));
        
        weights_additional_spectra.push_back(retspec);
        delete bbspec;

        return 0;
    }
    
    double interpol_fac = 0, diff;
    spectrum lowT_lowM, lowT_hiM, hiT_lowM, hiT_hiM;
    spectrum lowM, hiM;
    spectrum *result;

    //
    // Find surrounding 8 best spectra 
    // Closest in logG, Teff and metall
    //

    vector<spectrum*> level1(8);
    vector<spectrum*> level2(4);
    vector<spectrum*> level3(2);
    vector<spectrum*> level4(1);
    
    double thisweight;
    vector<double> interpol_weights(8);
    for (int i=0; i<interpol_weights.size(); i++) {
        interpol_weights.at(i) = 1;
    }

    vector<int> index_in_library;
    index_in_library.resize(8);
    for (int i=0; i<index_in_library.size(); i++) {
        index_in_library.at(i) = -1;
    }
    
    vector<double> interpol_weights_1(8);
    vector<double> interpol_weights_2(4);
    vector<double> interpol_weights_3(2);
    vector<double> interpol_weights_4(1);

    if (by_weight) {
        spectrum* tmp_spec;
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, -1, -1, -1, +1, &index_in_library.at(0));
        level1.at(0) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, -1, -1, +1, +1, &index_in_library.at(1));
        level1.at(1) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, -1, +1, -1, +1, &index_in_library.at(2));
        level1.at(2) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, -1, +1, +1, +1, &index_in_library.at(3));
        level1.at(3) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, +1, -1, -1, +1, &index_in_library.at(4));
        level1.at(4) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, +1, -1, +1, +1, &index_in_library.at(5));
        level1.at(5) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, +1, +1, -1, +1, &index_in_library.at(6));
        level1.at(6) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
        tmp_spec = find_spec(metall, logg, teff, alpha_enh, +1, +1, +1, +1, &index_in_library.at(7));
        level1.at(7) = (tmp_spec == NULL) ? 0 : new spectrum_blank(tmp_spec);
    } else {
        level1.at(0) = find_spec(metall, logg, teff, alpha_enh, -1, -1, -1, +1, &index_in_library.at(0));
        level1.at(1) = find_spec(metall, logg, teff, alpha_enh, -1, -1, +1, +1, &index_in_library.at(1));
        level1.at(2) = find_spec(metall, logg, teff, alpha_enh, -1, +1, -1, +1, &index_in_library.at(2));
        level1.at(3) = find_spec(metall, logg, teff, alpha_enh, -1, +1, +1, +1, &index_in_library.at(3));
        level1.at(4) = find_spec(metall, logg, teff, alpha_enh, +1, -1, -1, +1, &index_in_library.at(4));
        level1.at(5) = find_spec(metall, logg, teff, alpha_enh, +1, -1, +1, +1, &index_in_library.at(5));
        level1.at(6) = find_spec(metall, logg, teff, alpha_enh, +1, +1, -1, +1, &index_in_library.at(6));
        level1.at(7) = find_spec(metall, logg, teff, alpha_enh, +1, +1, +1, +1, &index_in_library.at(7));
    }

    /*
    // Find closest, assuming that metallicities already match, so we
    // can pick the second set of spectra
    //
    // Simulates in part the method of I_Isochrone
    // 
    double diff_Teff, diff_logg, totaldiff, mindiff = 100, min_Teff = 100, min_logg = 100;
    int best = -1;
    for (int i=4; i<8; i++) {
        if (level1.at(i) == NULL) continue;
        diff_Teff = abs(log10(teff) - log10(level1.at(i)->Teff));
        diff_logg = abs(logg - level1.at(i)->logg);
        totaldiff = diff_Teff; // + diff_logg;

        if (diff_Teff < 0.1 && diff_Teff < min_Teff) {
            min_Teff = diff_Teff;
            best = i;
            mindiff = totaldiff;
        }
    }
    for (int i=4; i<8; i++) {
        if (level1.at(i) == NULL) continue;
        diff_Teff = abs(log10(teff) - log10(level1.at(i)->Teff));
        diff_logg = abs(logg - level1.at(i)->logg);

        if (diff_Teff <= min_Teff && diff_logg < min_logg) {
            min_Teff = diff_Teff;
            min_logg = diff_logg;
            mindiff = totaldiff;
            best = i;
        }
    }
    if (best > 0) {
        double* pweight = weight(index_in_library.at(best));
        if (pweight != NULL) {
            *pweight += target_flux / level1.at(best)->get_normflux(normfilter)
                * additional_weight_factor;
        }
    } else {
        return 0;
    }
    */
    
    // cout << "Starting interpolating ..." << endl;
    
    // cout << "Searching for M=" << metall << " G=" << logg << " T=" << teff << endl;    
    //
    // for all 4 blocks
    // first interpolate in Teff
    //
    for (int i=0; i<4; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level2.at(i) = interpolate(level1.at(low), level1.at(hi), "Teff", teff, normfilter, &thisweight);
        for (int r=0; r<1; r++) {
            interpol_weights.at(low*1+r) *= (thisweight <= 0) ? 0 : thisweight;
            interpol_weights.at(hi*1+r) *= (thisweight <= 0) ? 0 : (1.0 - thisweight);
        }
        interpol_weights_1.at(low) = (thisweight < 0) ? 0 : 1.0 - thisweight;
        interpol_weights_1.at(hi)  = (thisweight < 0) ? 0 : thisweight;
        // cout << low << " " << hi << " --> " << interpol_weights_1.at(low) << " " << interpol_weights_1.at(hi) << endl;
    }
    // cout << endl;
    
    //
    // for all 2 blocks
    // next interpolate in logg
    // 
    for (int i=0; i<2; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level3.at(i) = interpolate(level2.at(low), level2.at(hi), "logg", logg, normfilter, &thisweight);
        for (int r=0; r<2; r++) {
            interpol_weights.at(low*2+r) *= (thisweight <= 0) ? 0 : thisweight;
            interpol_weights.at(hi*2+r) *= (thisweight <= 0) ? 0 : (1.0 - thisweight);
        }
        interpol_weights_2.at(low) = (thisweight < 0) ? 0 : 1.0 - thisweight;
        interpol_weights_2.at(hi)  = (thisweight < 0) ? 0 : thisweight;
        // cout << low << " " << hi << " --> " << interpol_weights_2.at(low) << " " << interpol_weights_2.at(hi) << endl;
    }
    // cout << endl;
    for (int i=0; i<4; i++) {
        if (level2.at(i) != 0) delete level2.at(i);
    }
    
    //
    // for all 2 blocks
    // next interpolate in metallicity
    // 
    for (int i=0; i<1; i++) {
        int low = 2*i + 0;
        int hi  = 2*i + 1;

        level4.at(i) = interpolate(level3.at(low), level3.at(hi), "metallicity", metall, normfilter, &thisweight);
        for (int r=0; r<4; r++) {
            interpol_weights.at(low*4+r) *= (thisweight <= 0) ? 0 : thisweight;
            interpol_weights.at(hi*4+r) *= (thisweight <= 0) ? 0 : (1.0 - thisweight);
        }
        interpol_weights_3.at(low) = (thisweight < 0) ? 0 : 1.0 - thisweight;
        interpol_weights_3.at(hi)  = (thisweight < 0) ? 0 : thisweight;
        // cout << low << " " << hi << " --> " << interpol_weights_3.at(low) << " " << interpol_weights_3.at(hi) << endl;
    }
    // cout << endl;
    for (int i=0; i<2; i++) {
        if (level3.at(i) != 0) delete level3.at(i);
    }

    if (level4.at(0) != 0) delete level4.at(0);

    double interpol_weights_all;
    for (int i=0; i<8; i++) {
        interpol_weights_all =
            interpol_weights_1.at(i / 1) * 
            interpol_weights_2.at(i / 2) *
            interpol_weights_3.at(i / 4);
/*XX        if (level1.at(i) == 0) {
            cout <<  "Spec_" << i << ": - weight=" << interpol_weights_all;
        } else {
            cout << "Spec_" << i
                 << ": M=" << fixed << setw(5) << setprecision(2) << level1.at(i)->metallicity
                 << " G=" << fixed << setw(5) << setprecision(2) << level1.at(i)->logg
                 << " T=" << fixed << setw(8) << setprecision(0) << level1.at(i)->Teff
                 << " weight=" << fixed << setprecision(4) << setw(8) << interpol_weights_all;
        }
        cout << " ("
             << fixed << setprecision(4) << setw(8) << interpol_weights_1.at(i / 1) << " "
             << fixed << setprecision(4) << setw(8) << interpol_weights_2.at(i / 2) << " "
             << fixed << setprecision(4) << setw(8) << interpol_weights_3.at(i / 4) << ")"
             << " ___" << level1.at(i)->filename << "___"
             << endl; 
             XX*/        
        if (index_in_library.at(i) >= 0) {
            double* pweight = weight(index_in_library.at(i));
            if (pweight != NULL) {
                *pweight += interpol_weights_all
                * target_flux / level1.at(i)->get_normflux(normfilter)
                * additional_weight_factor;
            }
        }
    }
  
    if (by_weight) {
        for (int i=0; i<level1.size(); i++) {
            if (level1.at(i) != 0) delete level1.at(i);
        }
    }

    return 0;
}


bool spec_library::is_valid()
{
    return well_prepared;
}


bool spec_library::test()
{
    double teff = 8500;
    double logg = 9.7;
    double metall = 1.8;

    spectrum *res = get_spec(metall, logg, teff, 0, "K", 0);

    cout << endl << endl << endl;
    cout << "   Requested: Teff=" << teff << ", logg=" << logg << ", metall=" << metall << endl;
    cout << "Interpolated: Teff=" << res->Teff << ", logg=" << res->logg
         << ", metall=" << res->metallicity << endl;

    return false;
}


bool spec_library::reset_weights()
{
    int bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2;
    vector< vector<teff_bin> > *select_m;
    vector< vector<logg_bin> > *select_t;

    // Reset all weights to zero
    for (bin_m = 0; bin_m < slibrary.size(); bin_m++) {
        for (bin_m2 = 0; bin_m2 < slibrary.at(bin_m).size(); bin_m2++) {
            select_m = &slibrary.at(bin_m).at(bin_m2).m;
            for (bin_t = 0; bin_t < select_m->size(); bin_t++) {
                for (bin_t2 = 0; bin_t2 < select_m->at(bin_t).size(); bin_t2++) {
                    select_t = &select_m->at(bin_t).at(bin_t2).t;
                    for (bin_g = 0; bin_g < select_t->size(); bin_g++) {
                        for(bin_g2 = 0; bin_g2 < select_t->at(bin_g).size(); bin_g2++) {
                            select_t->at(bin_g).at(bin_g2).weight = 0;
                        }
                    }                    
                }
            }
        }
    }

    // Also delete all additional spectra that might be around
    for (int i=0; i<weights_additional_spectra.size(); i++) {
        delete weights_additional_spectra.at(i);
    }
    weights_additional_spectra.clear();
    return true;
}

spectrum* spec_library::integrate_weights_to_spectrum()
{
    spectrum* fullspec = NULL;
    spectrum* addspec;
    spectrum* newspec;

    int bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2;
    vector< vector<teff_bin> > *select_m;
    vector< vector<logg_bin> > *select_t;
    double thisweight;

    // cout << "Integrating spectrum from weights ..." << endl;

    // Start with valid spectrum with all fluxes set to zero
    // This prevents the case of an illegal return spectrum in the case
    // that all weights are <= 0.
    fullspec = new spectrum(&wavelength, NULL);
    
    for (bin_m = 0; bin_m < slibrary.size(); bin_m++) {
        for (bin_m2 = 0; bin_m2 < slibrary.at(bin_m).size(); bin_m2++) {
            select_m = &slibrary.at(bin_m).at(bin_m2).m;
            for (bin_t = 0; bin_t < select_m->size(); bin_t++) {
                for (bin_t2 = 0; bin_t2 < select_m->at(bin_t).size(); bin_t2++) {
                    select_t = &select_m->at(bin_t).at(bin_t2).t;
                    for (bin_g = 0; bin_g < select_t->size(); bin_g++) {
                        for(bin_g2 = 0; bin_g2 < select_t->at(bin_g).size(); bin_g2++) {

                            // Skip spectra that do not contribute to integrated spec.
                            if (select_t->at(bin_g).at(bin_g2).weight <= 0) continue;
                            //cout << "Adding weight=" << select_t->at(bin_g).at(bin_g2).weight << endl;

                            /*XX
                            cout << fixed << setprecision(0) << setw(12) << encode_library_pos(bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2) << " "
                                 << fixed << setw(6) << setprecision(0) << select_t->at(bin_g).at(bin_g2).spec->Teff << " "
                                 << fixed << setw(6) << setprecision(2) << select_t->at(bin_g).at(bin_g2).spec->logg << " "
                                 << fixed << setw(6) << setprecision(1) << select_t->at(bin_g).at(bin_g2).spec->metallicity << " "
                                 << scientific << setprecision(5) << select_t->at(bin_g).at(bin_g2).weight << endl;XX*/

                            // Compute new, weighted contribution and add to stack
                            newspec
                                = *select_t->at(bin_g).at(bin_g2).spec
                                * select_t->at(bin_g).at(bin_g2).weight;
                            addspec = *fullspec + newspec;

                            // Delete the old intermediate spectra to avoid memory problems
                            // Future: Implement a += operator for spectra
                            delete fullspec;
                            delete newspec;
                            
                            fullspec = addspec;
                            // select_t->at(bin_g).at(bin_g2).spec->unload();
                        }
                    }                    
                }
            }
        }
    }

    for (int ms=0; ms<weights_additional_spectra.size(); ms++) {
        // similar to above, add weighted spectrum and then
        // delete the intermediate spectra
        addspec = *fullspec + weights_additional_spectra.at(ms);
        delete fullspec;
        fullspec = addspec;
    }
    
    return fullspec;
}


spectrum* spec_library::get_spec(int spec_id)
{
    return library.at(spec_id);
}

int spec_library::initialize_library_struct()
{

    slibrary.resize(4);
/*    for (int b_teff=0; b_teff<slibrary.size(); b_teff++) {
        slibrary.at(b_teff).resize(10);
        for (int b_logg=0; b_logg < slibrary.at(b_teff).size(); b_logg++) {
            slibrary.at(b_teff).at(b_logg).clear();
        }
    }
*/
    
    for (int cur_metall = 0; cur_metall < slibrary.size(); cur_metall++) {
        slibrary.at(cur_metall).clear();
    }
    
    return true;
}


int spec_library::encode_library_pos(int bin_m, int bin_m2, int bin_t,
                                     int bin_t2, int bin_g, int bin_g2)
{
    
    // Position is encoded as 32-bit number as follows:
    // Bit  1- 3: bin_m    (max. 8)
    // Bit  4- 8: bin_m2   (max. 32)
    // Bit  9-14: bin_t    (max: 64)
    // Bit 15-22: bin_t2   (max: 256)
    // Bit 23-26: bin_g    (max: 16)
    // Bit 27-32: bin_g2   (max: 64)
    // --> Maximum resolutions (assuming equal spacing):
    //     Metallicity: 0.125 dex
    //            Teff: 4 K
    //            logg: 0.0167 dex
    int position = 0 |
        ((bin_m & 0x07) << 29) | ((bin_m2 & 0x1f) << 24) |
        ((bin_t & 0x3f) << 18) | ((bin_t2 & 0xff) << 10) |
        ((bin_g & 0x0f) <<  6) | ((bin_g2 & 0x3f) <<  0);

    return position;
}



int spec_library::decode_library_pos(int pos, int* bin_m, int* bin_m2, int* bin_t,
                                     int* bin_t2, int* bin_g, int* bin_g2)
{

    if (bin_m  != NULL) *bin_m  = (pos >> 29) & 0x07;
    if (bin_m2 != NULL) *bin_m2 = (pos >> 24) & 0x1f;

    if (bin_t  != NULL) *bin_t  = (pos >> 18) & 0x3f;
    if (bin_t2 != NULL) *bin_t2 = (pos >> 10) & 0xff;

    if (bin_g  != NULL) *bin_g  = (pos >>  6) & 0x0f;
    if (bin_g2 != NULL) *bin_g2 = (pos >>  0) & 0x3f;

    return 0;
}


double* spec_library::weight(int index)
{
    int bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2;
    decode_library_pos(index,
                       &bin_m, &bin_m2, &bin_t, &bin_t2, &bin_g, &bin_g2);

    // Check if the coordinates are valid
    
    if (bin_m < 0 || bin_m2 < 0 || bin_t < 0 || bin_t2 < 0 || bin_g < 0 || bin_g2 < 0) return NULL;
    
    if (bin_m  >= slibrary.size()) return NULL;
    if (bin_m2 >= slibrary.at(bin_m).size()) return NULL;
    if (bin_t  >= slibrary.at(bin_m).at(bin_m2).m.size()) return NULL;
    if (bin_t2 >= slibrary.at(bin_m).at(bin_m2).m.at(bin_t).size()) return NULL;
    if (bin_g  >= slibrary.at(bin_m).at(bin_m2).m.at(bin_t).at(bin_t2).t.size()) return NULL;
    if (bin_g2 >= slibrary.at(bin_m).at(bin_m2).m.at(bin_t).at(bin_t2).t.at(bin_g).size()) return NULL;
    
    // All checks passed, so return the pointer to the weight
    return &slibrary.at(bin_m).at(bin_m2).m.at(bin_t).at(bin_t2).t.at(bin_g).at(bin_g2).weight;
    
}


bool spec_library::dump_structure()
{
    int bin_m, bin_m2, bin_t, bin_t2, bin_g, bin_g2;
    vector< vector<teff_bin> > *select_m;
    vector< vector<logg_bin> > *select_t;

    // Reset all weights to zero
    for (bin_m = 0; bin_m < slibrary.size(); bin_m++) {
        for (bin_m2 = 0; bin_m2 < slibrary.at(bin_m).size(); bin_m2++) {
            select_m = &slibrary.at(bin_m).at(bin_m2).m;
            for (bin_t = 0; bin_t < select_m->size(); bin_t++) {
                for (bin_t2 = 0; bin_t2 < select_m->at(bin_t).size(); bin_t2++) {
                    select_t = &select_m->at(bin_t).at(bin_t2).t;
                    for (bin_g = 0; bin_g < select_t->size(); bin_g++) {
                        for(bin_g2 = 0; bin_g2 < select_t->at(bin_g).size(); bin_g2++) {
                            cout << setw(2) << fixed << setprecision(0) << bin_m << " "
                                 << setw(2) << fixed << setprecision(0) << bin_m2 << " "
                                 << setw(2) << fixed << setprecision(0) << bin_t << " "
                                 << setw(2) << fixed << setprecision(0) << bin_t2 << " "
                                 << setw(2) << fixed << setprecision(0) << bin_g << " "
                                 << setw(2) << fixed << setprecision(0) << bin_g2 << "   :   "
                                 << setw(5) << fixed << setprecision(2) << select_t->at(bin_g).at(bin_g2).spec->metallicity << " "
                                 << setw(5) << fixed << setprecision(0) << select_t->at(bin_g).at(bin_g2).spec->Teff << " "
                                 << setw(5) << fixed << setprecision(2) << select_t->at(bin_g).at(bin_g2).spec->logg << endl;
                        }
                    }
                    cout << endl;
                }
            }
            cout << endl << endl;
        }
    }

    return true;
}


string spec_library::info()
{
    ostringstream oss;

    if (wavelength.size() > 0) {
        oss << "Spectral range: " << wavelength.at(0) << " - "
            << wavelength.at(wavelength.size()-1) << " covered by "
            << wavelength.size() << " datapoints";
    } else {
        oss << "Error: No wavelength datapoints read!";
    }
    
    return oss.str();
}
