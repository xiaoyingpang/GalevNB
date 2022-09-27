/* 
***********************************************************

  Spectra class to be used as part of the GALEV package
  adapted for Blackbody spectra

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECTRA_BLACKBODY_H__
#define __SPECTRA_BLACKBODY_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;

#include "filter.h"
#include "spectra.h"


class spectrum_blackbody : public spectrum
{

 protected:
    bool read_specfile();

 public:
    spectrum_blackbody(double temp, vector<double> *wl,
                       filter_library * filtlib, void* sl);

};


#endif
