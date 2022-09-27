

#ifndef __INDICES_FITFUNCT_H__
#define __INDICES_FITFUNCT_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include "indices.h"

using namespace std;
    
class indices_fitfunct : public indices
{
  protected:
    int compute_indices(spectrum* spec);
//    vector<index_data> indices;
    
  public:
    indices_fitfunct();
    indices_fitfunct(vector<bandpass>* bp_defs, spectrum* spec);

};



#endif
