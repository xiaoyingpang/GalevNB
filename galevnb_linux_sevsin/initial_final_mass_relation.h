


#ifndef __INITIAL_FINAL_MASS_RELATION_H__
#define __INITIAL_FINAL_MASS_RELATION_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "imf.h"

struct initial_final_mass_relation_entry
{
    double mass_min;
    double mass_max;
    double constant;
    double scaling;
    double offset;
};

class initial_final_mass_relation
{
  protected:
    vector<initial_final_mass_relation_entry> entry;
    string config_filename;
    
  public:
    initial_final_mass_relation(string configfile);
    ~initial_final_mass_relation();

    double get_mass_in_remnants(cimf* imf, double mass_min, double mass_max);

};



#endif
