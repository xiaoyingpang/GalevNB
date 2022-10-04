#include "filterlib.h"

filter_library::filter_library(string filename)
{

    return;
}

filter_library::filter_library()
{


    return;
}


int filter_library::add_filter(filter* ff)
{
    filters.push_back(ff);
    return filters.size();
}

int filter_library::add_filter(string filename, string name, double zpVega, double zpAB, double zpST)
{

    return filters.size();
}

int filter_library::add_filter(string filename, string name)
{
    filter* tmp = new filter(filename, name);
    return add_filter(tmp);
}

filter* filter_library::get_filter(int id)
{
    if (id >= filters.size()) return filters.at(0);
    
    return filters.at(id);
}

filter* filter_library::get_filter(string fname)
{
    for (int cf=0; cf<filters.size(); cf++) {
        if (filters.at(cf)->get_name() == fname) {
            return filters.at(cf);
        }
    }
    
    return NULL; //filters.at(0);
}

int filter_library::get_number_filters()
{
    return filters.size();
}
