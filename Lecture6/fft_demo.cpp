#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;

// Return Discrete Fourier Transform (DFT) of a complex data vector
vector< complex<double> >
discrete_transform(const vector< complex<double> >& data)
{
    int N = data.size();
    vector< complex<double> > transform(N, 0.0);
    const double pi = 4 * atan(1.0);
    const complex<double> i(0.0, 1.0);
    for (int k = 0; k < N; k++)
        for (int j = 0; j < N; j++) {
            double angle = 2 * pi * k * j / double(N);
            transform[k] += data[j] * exp(i * angle);
        }
    return transform;
}

// Return Fast Fourier Transform (FFT) using Danielson-Lanczos Lemma
vector< complex<double> >
recursive_transform(const vector< complex<double> >& data)
{
    int N = data.size();
    if (N == 1)             // transform is trivial
        return data;
    else if (N % 2 == 1)    // transform is odd, lemma does not apply
        return discrete_transform(data);

    // perform even-odd decomposition and transform recursively
    vector< complex<double> > even(N / 2, 0.0);
    for (int j = 0; j < N / 2; j++)
        even[j] = data[2 * j];
    even = recursive_transform(even);
    vector< complex<double> > odd(N / 2, 0.0);
    for (int j = 0; j < N / 2; j++)
        odd[j] = data[2 * j + 1];
    odd = recursive_transform(odd);
    const double pi = 4 * atan(1.0);
    const complex<double> i(0.0, 1.0);
    complex<double> W = exp(2.0 * i * pi / double(N));
    complex<double> Wk = 1.0;
    vector< complex<double> > transform(N, 0.0);
    for (int k = 0; k < N; k++) {
        transform[k] = even[k % (N / 2)] + Wk * odd[k % (N / 2)];
        Wk *= W;
    }
    return transform;
}

// Return Inverse Fast Fourier Transform (IFFT) using Danielson-Lanczos Lemma
vector< complex<double> >
recursive_inverse_transform(const vector< complex<double> >& data)
{
  std::vector<std::complex<double> > transform (data.size());
  std::copy( data.begin(), data.end(), transform.begin() );

  // conjugate the complex numbers
  for ( std::vector< std::complex<double> >::iterator itransform = transform.begin(),
	  itransformEnd = transform.end(); 
	itransform != itransformEnd; ++itransform ){
    *itransform = std::conj(*itransform);
  }

 
  // forward fft
  transform = recursive_transform( transform );
 
  // conjugate the complex numbers again
  for ( std::vector< std::complex<double> >::iterator itransform = transform.begin(),
	  itransformEnd = transform.end(); 
	itransform != itransformEnd; ++itransform ){
    *itransform = std::conj(*itransform);
  }

 
  // scale the numbers
  // conjugate the complex numbers
  for ( std::vector< std::complex<double> >::iterator itransform = transform.begin(),
	  itransformEnd = transform.end(); 
	itransform != itransformEnd; ++itransform ){
    *itransform /= transform.size();
  }

  return transform;
}

// Return one-sided power spectrum of transformed data
vector< vector<double> > power_spectrum(const vector< complex<double> >& data)
{
    int N = data.size();
    const double pi = 4 * atan(1.0);
    double omega = pi / N;
    vector< vector<double> > power(N / 2 + 1, vector<double>(2, 0.0));
    power[0][0] = 0.0;
    power[0][1] = pow(abs(data[0]), 2.0);
    for (int k = 1; k < N / 2; k++) {
        power[k][0] = k * omega;
        power[k][1] = pow(abs(data[k]), 2.0) + pow(abs(data[N - k]), 2.0);
    }
    if (N % 2 == 0) {
        power[N / 2][0] = (N / 2) * omega;
        power[N / 2][1] = pow(abs(data[N / 2]), 2.0);
    }
    for (int k = 0; k < power.size(); k++)
        power[k][1] /= pow(N, 2.0);
    return power;
}


int main()
{



    const double pi = 4 * atan(1.0);
    const complex<double> i(0.0, 1.0);
    const int N = 256;                   // number of transform points
 
    vector< complex<double> > data(N);

    string file_name("co2_mm_mlo.txt");
    // read the data file and
    ifstream data_file(file_name.c_str());
    if (data_file.fail()) {
        cerr << "cannot open " << file_name << endl;
	return 0;
    } else
      cout << " reading data file: " << file_name << endl;

    string line; 
    int j = 0; 
    while (getline(data_file, line) && j < N) {
        if (line.c_str()[0] != '#') {    
	  stringstream sline(line);
	  int year=0,month=0,days=-1;
	  double date=0., average = 0., interp = 0., trend = 0.;
	  sline >> year >> month >> date >> average >> interp >> trend >> days;
	  //std::cout << "read : date, interp = " << date << ", " << interp << endl;

	  data[j] = interp; 
	  ++j; 
	}
    }
    data_file.close();

    vector< complex<double> > P = recursive_transform(data);

    vector< complex<double> > Pinv = recursive_inverse_transform( P );

    ofstream file; 
    file_name = "fft_co2.data";
    file.open(file_name.c_str());
    for (int j = 0; j < P.size(); j++) {
        double ang_freq = j;
        double mag = abs(P[j]);
	double inv = abs(Pinv[j]);
        file << ang_freq << '\t' << data[j].real() << '\t' << mag << '\t' << inv << '\n';
    }
    file.close();
    cout << " wrote P and Pinv values in " << file_name << endl;
}
