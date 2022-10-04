/* 
***********************************************************

  Filter class to be used as part of the GALEV package

  written by: Ralf Kotulla
  (c) 2008, The GALEV Team

***********************************************************
*/

#ifndef __FILTER_H__
#define __FILTER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


using namespace std;


struct datapoint
{
    double wavelength;
    double response;
};

class filter
{
 protected:
//    vector<double> wavelength;
//    vector<double> response;
    double minwavel, maxwavel;
    int nwavel;

    std::vector<datapoint> data;

    vector<double> optimized_response;
    vector<double> resolution_element_width;
    bool optimized;
    double integrated_width;
    
    double zpv; 
    double zpab; 
    double zpst;

    string name;
    string filename;
    
 public:
    filter();
    filter(string filename, string name);
    ~filter();
    
    bool optimize(vector<double>* lambda);
    
    double get_response(double wavelength);
    int dump();
    double lambda_start();
    double lambda_end();
    string get_name();

    double get_response_curve(vector<double>** responses,
                              vector<double>** element_width);
    
};

#endif
