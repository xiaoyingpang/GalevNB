/* 
***********************************************************

  IMF class to be used as part of the GALEV package

  This class can deal with binned IMFs. The number of bins
  is not limited to allow maximum flexibility.

  Necessary parameters (obtained from config class):
  - IMF_FILE
  - IMF_NORMRANGE
  
  written by: Peter Anders
  (c) 2011, The GALEV Team

***********************************************************
*/


#ifndef __IMF_BINNED_H__
#define __IMF_BINNED_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "imf.h"
#include "config.h"

struct binned_part
{
    double lower_mass_limit;
    double upper_mass_limit;
    double star_count;
};

class cimf_binned : public cimf
{
  private:
//    double integrate_binned(double min, double max);

  protected:
    double normalize_minmass, normalize_maxmass;
    vector<binned_part> IMF_Bins;
    
    
  public:
    cimf_binned(cconfig* cc);
    ~cimf_binned();
    
    double get_maxmass();
    double get_minmass();
    double integrate_imf(double ml, double mu, double mass_exponent);
};

#endif
