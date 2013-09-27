#ifndef CPT_STD_HPP
#define CPT_STD_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

namespace cpt {

static void cpt_error(const char* s)
{
    cerr << " cpt_error: " << s << endl;
    exit(EXIT_FAILURE);
}

static void cpt_error(string s) { cpt_error(s.c_str()); }

}  // end namspace cpt

#endif /* CPT_STD_HPP */
