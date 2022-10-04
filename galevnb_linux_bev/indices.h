

#ifndef __INDICES_H__
#define __INDICES_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include "spectra.h"

using namespace std;

struct bandpass
{
    string name;
    int    mode;
    
    double cont_ll;
    double cont_lr;
    double pass_l;
    double pass_r;
    double cont_rl;
    double cont_rr;

    vector<double> add_params;
};

struct index_data
{
    double cont_l;
    double flux;
    double cont_r;
};

    
class indices
{
  protected:
    vector<bandpass>* bandpass_defs;
    vector<index_data> indices_data;
    virtual int compute_indices(spectrum* spec);
    
    
  public:
    indices();
    indices(vector<bandpass>* bp_defs, spectrum* spec);
    indices(vector<bandpass>* bp_defs);
    
    indices* operator+(indices* index);
    indices* operator-(indices* index);
    indices* operator*(double factor);
    indices* operator/(double factor);
};



#endif
