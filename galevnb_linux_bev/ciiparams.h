
#ifndef __CIIPARAMS_H__
#define __CIIPARAMS_H__

#include "spectra.h"
#include "speclib.h"
#include "isochrone.h"
#include "imf.h"
#include "lycphotons.h"
#include "emissionlines.h"
#include "config.h"

struct ciiparams
{
    spectrum** ret_isospec;

    spec_library* speclib;
    isochrone* isos;
    cimf* imf;
    lycphotons* lyc;
    emissionlines* emlines;
    
    int metall;
    int age;

    cconfig* cconf;

    bool first_iso;
    bool last_iso;
};

#endif
