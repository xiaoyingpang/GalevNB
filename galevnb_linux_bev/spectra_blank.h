/* 
***********************************************************

  Spectra class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECTRA_BLANK_H__
#define __SPECTRA_BLANK_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;

#include "filterlib.h"
#include "filter.h"
#include "spectra.h"

class spectrum_blank : public spectrum
{

 public:
    spectrum_blank();

    spectrum_blank(spectrum *spec);
    ~spectrum_blank();

};

#endif
