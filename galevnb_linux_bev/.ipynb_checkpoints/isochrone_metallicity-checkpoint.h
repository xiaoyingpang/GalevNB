
#ifndef __ISOCHRONE_METALLICITY_H__
#define __ISOCHRONE_METALLICITY_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "constants.h"
#include "isochrone_age.h"

class isochrone_metallicity
{
  private:
    vector<string> filenames;
    vector<int> isoIDs;
    vector<bool> file_loaded;
    
    int read_isofile(string filename, int isoid, int metall_id, double metallicity, double alpha_enh);
    bool autoload;
    
    
  public:
    isochrone_metallicity(bool p_autoload);
    isochrone_metallicity();
    ~isochrone_metallicity();
    int add_file(string filename, int isoID);
    int load();
    int unload();
    int search_isochrone_age(double age);

    bool loaded;
    double metallicity;
    double alpha_enh;
    
    int metall_id;
    vector<isochrone_age> ages;    
};


#endif
