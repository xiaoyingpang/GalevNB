/* 
***********************************************************

  Spectra class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __SPECTRA_BINARY_H__
#define __SPECTRA_BINARY_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;

#include "spectra.h"
#include <stdint.h>

typedef int32_t int32;


class spectrum_binary : public spectrum
{
 protected:
    virtual bool read_specfile();
    virtual bool write_specfile(string filename);
    
 public:
    spectrum_binary();
    spectrum_binary(spectrum* spec);
    ~spectrum_binary();
};

struct SPECTRUM_BINARY_FILEHEAD 
{
    char    id_label[12];           // GALEVSPEC
    int32 version;               // Version to allow future upgrades
    int32 number_blocks;         // default: 1
    int32 startpos;              //
    char    reserved[76];           // to make this block 100 bytes long for future updates
};

struct SPECTRUM_BINARY_BLOCKHEAD
{
    int32 head_length;           // for now this should be 200
    int32 block_type;            // spectrum (from library, isoint, or GALEV), wavelength definition, etc
    int32 startpos_next_block;   // where in the file does the next block start
    char reserved[188];          // fill up the 200 bytes to have space for future updates
};
enum SPECTRUM_BINARY_BLOCKTYPES
{
    SPECBINTYPE_SPECLIB = 1,
    SPECBINTYPE_ISOSPEC,
    SPECBINTYPE_GALEVSPEC
};

struct SPECTRUM_BINARY_SPECLIB
{
    int32_t number_lambdas;        // how many lambdas in this spec
    int32_t number_calibfilters;
};
struct SPECTRUM_BINARY_SPECLIB_CALIBFILTER
{
    char name[12];
    double magnitude;
};


#endif
