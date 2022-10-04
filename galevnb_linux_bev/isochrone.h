/* 
***********************************************************

  Isochrone class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __ISOCHRONE_H__
#define __ISOCHRONE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

using namespace std;

#include "constants.h"
#include "imf.h"

// Includes for derived isochrone classes have to be
// inserted at the end of this file, otherwise results in mutual includes
// that can not be compiled

#define ISOID_PADOVA    1
#define ISOID_GENEVA    2
#define ISOID_BINARIES  3
#define ISOID_PADOVA_WEBCMD_V23  4

#include "isochrone_metallicity.h"
#include "isochrone_age_mass.h"
#include "isochrone_padova.h"
#include "isochrone_padova_webcmd_v23.h"
#include "isochrone_binaries.h"

class isochrone
/**
 *
 * @brief Isochrone class
 *
 * Handle everything concerning the isochrones as a whole, i.e.
 * to make sure the right files are read and the different metallicities
 * are combined into one set.
 *
 * */
{

 private:
  
 public:
    //! Store pointers to all metallicity classes
    vector<isochrone_metallicity*> metalls;

    isochrone();
    ~isochrone();
    int search_isochrone_age(int metall, double age);
    int read_isochrone_set(string listfile, string isoname);
    isochrone_age_mass* find(int metall_id, int isoage, int massbin);

};



#endif
