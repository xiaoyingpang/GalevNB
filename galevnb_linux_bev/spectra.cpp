#include "spectra.h"
#include "speclib.h"

spectrum::spectrum()
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


spectrum::~spectrum()
{

    //    delete &flux;
    //    delete &norm_fluxes;
    unload();
    
    return;
}

void spectrum::destroy()
{
//    this->~spectrum();
//    cout << "Trying to delete spectrum at " << this << endl;
    delete this;
//    spec_counter--;
//    cout << spec_counter << " ";
//    norm_fluxes.clear();
//    flux.clear();     
}


bool spectrum::read_specfile()
{
//  cout << "Trying to access base class (" << loaded << "==" << true << ")" << endl;
    return loaded;
}

spectrum::spectrum(string p_filename, 
         double p_teff, double p_logg, double p_metall, double p_alpha_enh,
         int p_skip_front, int p_read, int p_skip_end, 
         vector<double>* p_wl, bool load_now,
         filter_library* filtlib, void* sl)
{
//    cout << "Creating base class spectrum #1" << endl;
    flux = NULL;
    unload();
    
    norm_fluxes.clear();
    norm_fluxes_names.clear();

    for (int i=0; i<norm_fluxes.size(); i++) {
        set_normflux(filterlib->get_filter(i)->get_name(), 0.0);
    }
    set_normflux("MBol", 0.0);
    
    spectrum(p_wl,sl);
    return;
}

spectrum::spectrum(vector<double>* p_wl, void* sl)
{
//    cout << "Creating base class spectrum #3" << endl;
    flux = NULL;
    unload();
    
    wavelength = p_wl;
    flux = new vector<double>;
    flux->resize(wavelength->size());

    Teff = -99;
    logg = -99;
    metallicity = -99;
    alpha_enh = -99;
    valid = true;

    filename = "";
    loaded = true;

    for (int i=0; i<flux->size(); i++) {
        flux->at(i) = 0.00;
    }
    // bolometric_flux = 0;

    return;
}


void spectrum::dump(int steps)
{
    if (steps <= 0) steps = 5;
    
    cout << "# Trying to dump file ..." << endl;
    if (read_specfile()) {
        cout << "# Dumping " << flux->size() << " entries ..." << endl;
        for (int i=0; i<flux->size(); i+=steps) {
            cout << setw(6) << fixed << setprecision(0) << i+1 << "  " 
                 << setw(12) << scientific << setprecision(4) << wavelength->at(i) << "  " 
                 << setw(12) << scientific << setprecision(4) << flux->at(i) << endl;
        }
        cout << "# Dumping done!" << endl;
    }
    return;
}

spectrum* spectrum::operator+(spectrum* spec)
{
    return add(spec);
};

spectrum* spectrum::operator-(spectrum* spec)
{
    spectrum* dum = spec->mult(-1.0);
    spectrum* out = add(dum);
    delete dum;
    
    return out;
};


spectrum* spectrum::operator*(double scalefac)
{
//    read_specfile();
//    spectrum *out; // = *this;

    return mult(scalefac);
};


spectrum* spectrum::mult(double scalefac)
{
    read_specfile();
    spectrum *out = new spectrum(wavelength, speclib);
//    spectrum *out = ((spec_library*)speclib)->new_spec();
      //new spectrum; // = *this;

    out->filename = "";
    out->loaded = true;
    out->wavelength = wavelength;
    out->valid = true;

    out->Teff = Teff;
    out->logg = logg;
    out->metallicity = metallicity;
    out->alpha_enh = alpha_enh;

    out->flux->resize(flux->size());
    for (int i=0; i<flux->size(); i++) {
        out->flux->at(i) = flux->at(i) * scalefac;
    }

    out->filterlib = filterlib;
    out->norm_fluxes.resize(norm_fluxes.size());
    out->norm_fluxes_names.resize(norm_fluxes_names.size());
    for (int i=0; i < out->norm_fluxes.size(); i++) {
        out->norm_fluxes.at(i) = norm_fluxes.at(i) * scalefac;
        out->norm_fluxes_names.at(i) = norm_fluxes_names.at(i);
    }

    return out; 
}


spectrum* spectrum::operator/(double scalefac)
{
    return mult(1/scalefac); 
};




spectrum* spectrum::add(spectrum* spec)
{
    bool read_spec_1, read_spec_2;
    
    read_spec_1 = this->read_specfile();
    read_spec_2 = spec->read_specfile();
    
    if (!read_spec_1 || !read_spec_2 ) {
        cout << "Error while adding spectra, return NULL pointer" << endl << flush;
        return NULL;
    }

    spectrum *out = new spectrum(wavelength, speclib);
    out->filename = "";
    out->loaded = true;
    out->wavelength = wavelength;
    out->valid = true;

    out->Teff = out->logg = out->metallicity = out->alpha_enh = -99;

    out->flux->resize(wavelength->size());
    for (int i=0; i<wavelength->size(); i++) {
        out->flux->at(i) = flux->at(i) + spec->flux->at(i);
    }

    out->norm_fluxes.clear();
    out->norm_fluxes_names.clear();
    for (int i=0; i<norm_fluxes.size(); i++) {
        string filter = norm_fluxes_names.at(i);
        double flux1 = spec->get_normflux(filter);
        double flux2 = get_normflux(filter);

        if (flux1 >= 0 && flux2 >= 0) {
            out->norm_fluxes.push_back(flux1+flux2);
            out->norm_fluxes_names.push_back(filter);
        }
    }

    return out; 
}


void spectrum::operator*=(const double factor)
{
    read_specfile();
    for (int i=0; i<flux->size(); i++) {
        flux->at(i) *= factor;
    }

    for (int i=0; i < norm_fluxes.size(); i++) {
        norm_fluxes.at(i) *= factor;
    }
    
    return;
}


void spectrum::operator=(spectrum spec)
{
    spec.read_specfile();
    
    loaded = true;
    valid = spec.valid;

    flux->resize(spec.extract_lambdas()->size());
    for (int i=0; i<flux->size(); i++) {
        flux->at(i) = spec.extract_fluxes()->at(i);
    }
    filename = spec.filename;
    wavelength = spec.wavelength;

    Teff = spec.Teff;
    metallicity = spec.metallicity;
    logg = spec.logg;
    alpha_enh = spec.alpha_enh;

    filterlib = spec.filterlib;
    norm_fluxes.resize(spec.norm_fluxes.size());
    norm_fluxes_names.resize(spec.norm_fluxes_names.size());
    for (int i=0; i < norm_fluxes.size(); i++) {
        norm_fluxes.at(i) = spec.norm_fluxes.at(i);
        norm_fluxes_names.at(i) = spec.norm_fluxes_names.at(i);
    }
    
    return;
}

double spectrum::operator*(filter ff)
{
    return integrate_flux(&ff);
}

double spectrum::integrate_flux(filter *ff)
{
    int min, max;
    double iflux, intf;
    read_specfile();

    iflux=0.0;
    intf=0.0;
    if((wavelength->at(0) > ff->lambda_end()) ||
       (wavelength->at(wavelength->size()-1) < ff->lambda_start())) {
        return (0);
    } else {

        vector<double>* response;
        vector<double>* width;

        ff->optimize(wavelength);
        intf = ff->get_response_curve(&response, &width);

        for (int i=0; i<wavelength->size(); i++) {
            // cout << ".";
            iflux += flux->at(i) * response->at(i) * width->at(i);
        }
        // cout << endl;
        
    }

    if(iflux>1e-99 && intf>0) {
        return(iflux/intf);
    } else {
        return(0);
    }

}

double spectrum::integrate_flux()
{
    double totalflux = 0, totalwidth = 0, width;
    
    int lower, upper;
    for (int i=0; i<wavelength->size(); i++) {
        lower = (i > 0) ? i-1 : 0;
        upper = (i < wavelength->size()-1) ? i+1 : wavelength->size() - 1;

        width = 0.5 * (wavelength->at(upper) - wavelength->at(lower));
        totalwidth += width;
        totalflux += flux->at(i) * width;
    }
    return (totalflux / totalwidth);
}

double spectrum::integrate_flux(double start, double end)
{
    if (start < 0) start = wavelength->at(0);
    if (end < 0) end = wavelength->at(wavelength->size()-1);
    
    double
        flux_total = 0,
        width = (end - start),
        flux_start, flux_end;
    
    for(int curl=0; curl<wavelength->size();curl++) {
        if (wavelength->at(curl+1) < start) {
            continue;
        } else if (wavelength->at(curl) > end) {
            break;
        } else {
            // wavelength(curl) < start && wavelength(curl+1) > start
            if (wavelength->at(curl+1) > end) {
                // both start and end within [curl;curl+1]
                flux_start = flux->at(curl) +
                    (flux->at(curl+1) - flux->at(curl)) *
                    ((start - wavelength->at(curl)) / (wavelength->at(curl+1) - wavelength->at(curl)));
                flux_end = flux->at(curl) +
                    (flux->at(curl+1) - flux->at(curl)) *
                    ((end - wavelength->at(curl)) / (wavelength->at(curl+1) - wavelength->at(curl)));
                flux_total += ((flux_start + flux_end) / 2) * (end - start);
            } else {
                flux_start = flux->at(curl) +
                    (flux->at(curl+1) - flux->at(curl)) *
                    ((start - wavelength->at(curl)) / (wavelength->at(curl+1) - wavelength->at(curl)));
                flux_end = flux->at(curl+1);
                flux_total += ((flux_start + flux_end) / 2) * (wavelength->at(curl+1) - start);
                start = wavelength->at(curl+1);
            }
        }
    }
    
    return (flux_total / width);
}


spectrum::spectrum(string p_filename, vector<double> *wl,
                   int *p_skip_begin, int *p_read, int *p_skip_end,
                   filter_library* filtlib, void* sl)
{

}


bool spectrum::is_valid()
{
    return valid;
}



vector<double>* spectrum::extract_lambdas() 
{
//    read_specfile();
    return wavelength;
}

vector<double>* spectrum::extract_fluxes()
{
    read_specfile();
    return flux;
}

double spectrum::get_normflux(string filter)
{
    read_specfile();
    for (int cf=0; cf<norm_fluxes.size(); cf++){
        if (norm_fluxes_names.at(cf) == filter) {
            return norm_fluxes.at(cf);
        }
    }
    return -1;
}


double spectrum::get_normflux(int filterid)
{
    cout << "using non-std. function" << endl;
    
    read_specfile();
    if (filterid < 0) {
        return (bolometric_flux <= 0 ? 1e-50 : bolometric_flux);
    } else if (filterid >= norm_fluxes.size()) {
        return -99;
    } 
    return (norm_fluxes.at(filterid) <= 0 ? 1e-50 : norm_fluxes.at(filterid));
}

    
int spectrum::set_normflux(string filter, double value)
{
    for (int cf=0; cf<norm_fluxes.size(); cf++){
        if (norm_fluxes_names.at(cf) == filter) {
            norm_fluxes.at(cf) = value;
            return cf;
        }
    }

    // If we happen to come here then this filter does not have
    // a normflux entry yet, so let's create a new one.
    norm_fluxes_names.push_back(filter);
    norm_fluxes.resize(norm_fluxes_names.size());
    norm_fluxes.at(norm_fluxes_names.size()-1) = value;
        
    return (norm_fluxes.size() - 1);
}


//
// Unload routine to free up memory
//
bool spectrum::unload()
{
    // Reset the fluxes, but leave the normfluxes untouched
    if (flux != NULL) {
        delete flux;
        flux = NULL;
    }
    
    loaded = false;
    return true;
}



//
// Unload routine to free up memory
//
string spectrum::info()
{

    cout << "Flux length: " << flux->size() << endl;
    cout << "Wavelength length: " << wavelength->size() << endl;
    
    return "";
}
