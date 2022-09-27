

/* 
***********************************************************

  filter library class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __FILTERLIB_H__
#define __FILTERLIB_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;

#include "filter.h"

class filter_library
{
  private:
    vector<filter*> filters;
    
  public:
    filter_library(string filename);
    filter_library();

    int add_filter(filter* ff);
    int add_filter(string filename, string name);
    int add_filter(string filename, string name, double zpVega, double zpAB, double zpST);

    filter* get_filter(int id);
    filter* get_filter(string fname);

    int get_number_filters();
    
};


    
    
#endif
