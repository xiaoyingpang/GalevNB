
#include "filter.h"

filter::filter()
{
    datapoint dp;
    nwavel = 40;
    
    double r[] = {
        0.000,0.006,0.030,0.060,0.134,0.302,0.567,0.841,0.959,
        0.983,0.996,1.000,0.996,0.987,0.974,0.957,0.931,0.897,
        0.849,0.800,0.748,0.698,0.648,0.597,0.545,0.497,0.447,
        0.397,0.345,0.297,0.252,0.207,0.166,0.129,0.095,0.069,
        0.043,0.024,0.009,0.000
    };

    double startwl=3600;
    for(int i=0; i<40 ; i++) {
        dp.wavelength = startwl + (i*50);
        dp.response = r[i];
        data.push_back(dp);
    }

    minwavel=data.at(0).wavelength;
    maxwavel=data.at(data.size()-1).wavelength;

    optimized = false;
}

filter::~filter()
{

}


filter::filter(string file, string fname)
{
    filename = file;
    
    // cout << "Reading filter-file " << filename << " ..." << endl << flush;
    ifstream filterfile(filename.c_str(), ios_base::in);
    if (!filterfile.good()) {
        cout << "Unable to open filter-file " << filename 
             << ", aborting ..." << endl;
    }
    string line;
    datapoint dp;

    minwavel = 1e99;
    maxwavel = -1e99;
    while(getline(filterfile, line, '\n')) {
        if (line.length() <= 0) continue;
        if (line.at(0) == '#') continue;
        
        istringstream iss(line);
        iss >> dp.wavelength >> dp.response;
        if (dp.wavelength >= maxwavel) {
            maxwavel = dp.wavelength;
            data.push_back(dp);
        } else {
            minwavel = dp.wavelength;
            data.insert(data.begin(), dp);
        }
    }

    dp.wavelength = data.at(0).wavelength - 1;
    minwavel = dp.wavelength;
    dp.response = 0;
    data.insert(data.begin(), dp);

    dp.wavelength = data.at(data.size()-1).wavelength + 1;
    maxwavel = dp.wavelength;
    dp.response = 0;
    data.push_back(dp);

    nwavel = data.size();

    name = fname;
    optimized = false;
}


double filter::get_response(double wavelength)
{
    int m;

    if(wavelength< minwavel) {
        return(0.0);
    }

    if(wavelength> maxwavel){
        return(0.0);
    }

    int li=0;
    int re=nwavel-1;
    m = ((li+re)/2);
    while(1) {
        if(wavelength >= data.at(m).wavelength && wavelength<=data.at(m+1).wavelength) 
            return(data.at(m).response+((data.at(m+1).response-data.at(m).response)
                                        /(data.at(m+1).wavelength-data.at(m).wavelength))*(wavelength-data.at(m).wavelength)); 
        if(wavelength < data.at(m).wavelength){
            re=m;
        }
        if(wavelength > data.at(m).wavelength) {
            li=m;
        }
        m=((li+re)/2);
    }

    return 0;
}


int filter::dump()
{
    cout << endl << endl << endl << endl
         << "# " << name << " (" << filename << ")" << endl;

    int i;
    for (i=0; i<data.size(); i++) {
        cout << fixed << setprecision(4) << data.at(i).wavelength << "  "
             << fixed << setprecision(5) << data.at(i).response << endl;
    }

    return i;
}


double filter::lambda_start()
{
    return minwavel;
}

double filter::lambda_end()
{
    return maxwavel;
}

string filter::get_name()
{
    return name;
}

bool filter::optimize(vector<double>* lambda)
{
    if (optimized) return true;

    integrated_width = 0;
    optimized_response.resize(lambda->size());
    resolution_element_width.resize(lambda->size());

    int lower, upper;
    for (int i=0; i<lambda->size(); i++) {

        optimized_response.at(i) = get_response(lambda->at(i));

        lower = (i > 0) ? i-1 : 0;
        upper = (i < lambda->size()-1) ? i+1 : lambda->size()-1;
        resolution_element_width.at(i) = 0.5 * (lambda->at(upper) - lambda->at(lower));

        integrated_width += optimized_response.at(i) * resolution_element_width.at(i);

        // cout << lambda->at(i) << " " << optimized_response.at(i) << " "
        //      << resolution_element_width.at(i) << " " << integrated_width << endl;
        
    }
    optimized = true;
    
    return optimized;
}

double filter::get_response_curve(vector<double>** responses,
                                  vector<double>** element_width)
{
    *responses = &optimized_response;
    *element_width = &resolution_element_width;
    return integrated_width;
}

