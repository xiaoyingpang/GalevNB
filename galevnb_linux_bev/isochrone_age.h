
#ifndef __ISOCHRONE_AGE_H__
#define __ISOCHRONE_AGE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "constants.h"
#include "isochrone_age_mass.h"

class isochrone_age
{

// private:

 public:
    isochrone_age();
    ~isochrone_age();

    int add_mass(isochrone_age_mass *iam);
    double age;

    double get_weight(cimf* imf, int index, int weightmode);
    vector<isochrone_age_mass*> masses;
};



#endif
