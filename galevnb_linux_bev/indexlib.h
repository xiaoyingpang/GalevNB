

#ifndef __INDEXLIB_H__
#define __INDEXLIB_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

#include "indices.h"
#include "indices_fitfunct.h"

using namespace std;

class index_library
{
  private:
    string definition_file;
    string bandpass_definitions;
    string index_mode;
    
    bool prepare_library();
    bool already_prepared;
    
    vector<bandpass> index_defs;
        
  public:
    index_library();
    index_library(string definition_file, string indexmode);
    bool getlist();
    indices* get_indices(spectrum* spec);

};



#endif
