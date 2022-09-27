

#ifndef __ISOCHRONE_AGE_MASS_H__
#define __ISOCHRONE_AGE_MASS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#define WEIGHTMODE_MASS   1
#define WEIGHTMODE_NSTARS 2

#include "constants.h"
#include "imf.h"

class isochrone_age_mass
{
 public:
    isochrone_age_mass() {};
    isochrone_age_mass(string line);
    ~isochrone_age_mass();
    
    virtual double get_mass_now() = 0;
    virtual double get_mass_zams() = 0;
    virtual double get_age() = 0;
    virtual double get_logage() = 0;

    virtual double get_Teff() = 0;
    virtual double get_logg() = 0;
    virtual double get_logL() = 0;
    virtual double get_magnitude(string filtername) = 0;

    virtual double get_weight(cimf* imf, int weightmode,
                              double prevmass, double nextmass);
};


#endif
