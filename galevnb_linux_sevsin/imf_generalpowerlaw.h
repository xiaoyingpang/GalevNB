/* 
***********************************************************

  IMF class to be used as part of the GALEV package

  This class can deal with broken powerlaws. The number of segments
  is not limited to allow maximum flexibility.

  Necessary parameters (obtained from config class):
  - IMF_PARAMS
  - IMF_NORMRANGE
  
  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/


#ifndef __IMF_GENERALPOWERLAW_H__
#define __IMF_GENERALPOWERLAW_H__

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

struct generalpowerlaw_part
{
    double lower_limit;
    double upper_limit;
    double slope;
    double scalefac;
};

class cimf_generalpowerlaw : public cimf
{
  private:
    bool match_parts();
    double integrate_generalpowerlaw(double min, double max, double deltaslope);

  protected:
    double normalize_minmass, normalize_maxmass;
    vector<generalpowerlaw_part> gpl_parts;
    
    
  public:
    cimf_generalpowerlaw(cconfig* cc);
    ~cimf_generalpowerlaw();
    
    double get_maxmass();
    double get_minmass();
    double integrate_imf(double ml, double mu, double delta_slope);
};

#endif
