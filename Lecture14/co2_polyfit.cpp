#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <polyfit.hpp>

int main()
{
  using namespace std; 


    // data downloaded from http://neic.usgs.gov/neis/epic/epic_global.html
    const string co2_data("co2_mm_mlo.txt");

    cout << " CO2 data " << endl;

    vector<double> x_values, y_values;

    // read the data file
    ifstream data_file(co2_data.c_str());
    if (data_file.fail()) {
        cerr << "cannot open " << co2_data << endl;
        exit(EXIT_FAILURE);
    } else
      cout << " reading data file: " << co2_data << endl;
    string line; 
    int j = 0; 
    while (getline(data_file, line)) {
        if (line.c_str()[0] != '#') {    
          stringstream sline(line);
          int year=0,month=0,days=-1;
          double date=0., average = 0., interp = 0., trend = 0.;
          sline >> year >> month >> date >> average >> interp >> trend >> days;
	  x_values.push_back( date );
	  y_values.push_back( interp );
          ++j; 
        }
    }

    data_file.close();




    cpt::Matrix<double,1> xvec(x_values.size());
    cpt::Matrix<double,1> yvec(x_values.size());
    cpt::Matrix<double,1> sigma(x_values.size());


    for ( unsigned int i = 0; i < x_values.size(); i++ ) {
      xvec[i] = x_values[i];
      yvec[i] = y_values[i];
      sigma[i] = 0.001; //sqrt(y_values[i]);
    }

    // cout << "Data are : " << endl;
    // cout << xvec << endl;
    // cout << yvec << endl;    
    // cout << sigma << endl;


    const int POLYNOMIAL_ORDER = 3;
    // fit histogram to straight line
    cpt::Matrix<double,1> a_fit(POLYNOMIAL_ORDER);
    cpt::Matrix<double,1> sig_a(POLYNOMIAL_ORDER);
    cpt::Matrix<double,1> yy(x_values.size());
    double chi2 = 0.0;
    polyfit(xvec, yvec, sigma, POLYNOMIAL_ORDER, a_fit, sig_a, yy, chi2);
    cout << "     a = " << a_fit << endl;
    cout << " sig_a = " << sig_a << endl;
    cout << "  chi2 = " << chi2 << endl;

    
    ofstream co2_polyfit_output("co2_polyfit_output.txt");
    char buff[1000];
    for ( unsigned int i = 0 ; i < x_values.size(); ++i ) {
      sprintf(buff, "%6.2f %6.2f %6.2f", xvec[i], yvec[i], yy[i]);
      co2_polyfit_output << buff << endl;
    }

    co2_polyfit_output.close();
}
