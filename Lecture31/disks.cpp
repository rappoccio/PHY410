#include "cptstd.hpp"
#include "graphics.hpp"
#include "linalg.hpp"
#include "random.hpp"
using namespace cpt;


class Disks {

public : 
  Disks(                    // set d_0, alpha, and position the disks
    double nu,              // auxiliary parameter varied from 0 to 7
    int rows=16,            // number of rows
    int columns=14) :       // number of columns
    A(1.0), N(224), r(N,2)
  {
    double d = 1 / double(columns);
    d_0 = d * (1 - pow(2.0, nu - 8));
    alpha = d - d_0;
    double dx = d;
    double dy = d * sqrt(3.0) / 2;
    for (int j = 0; j < rows; j++) {
        for (int i = 0; i < columns; i++) {
            int n = i + j * columns;
            r[n][0] = i * dx + (j % 2) * dx / 2.0;
            r[n][1] = j * dy;
       }
    }

  }

  bool overlap(                           // return true if overlap, else false
	       const Matrix<double,1>& r1,         // position vector of first disk
	       const Matrix<double,1>& r2)         // position vector of second disk
  {
    Matrix<double,1> dr = r1 - r2;
    // if separation in x or y is > A/2 check overlap with closest image
    for (int k = 0; k < 2; k++)
      if (abs(dr[k]) > 0.5 * A)
	dr[k] *= 1 - 1 / abs(dr[k]);
    return dot(dr, dr) < d_0 * d_0;
  }

  bool metropolis_step(int i) // apply Metropolis algorithm to disk i
  {
    Matrix<double,1> r_trial = r[i];
    for (int k = 0; k < 2; k++) {
      r_trial[k] += alpha * (2 * rng.rand() - 1.0);
      // apply periodic boundary conditions
      if (r_trial[k] < 0)
	r_trial[k] += A;
      if (r_trial[k] >= A)
	r_trial[k] -= A;
    }
    // accept move if there is no overlap with any other disk
    bool accept = true;
    for (int j = 0; j < N; j++)
      if (j != i && overlap(r_trial, r[j])) {
	accept = false;
	break;
      }
    if (accept)
      r[i] = r_trial;
    return accept;
  }

  void monte_carlo_step()
  {
    for (int i = 0; i < N; i++)
      metropolis_step(i);
  }


  void radial_distribution_histogram(
				     double nu,
				     double K = 1.5,
				     int number_of_bins = 64,
				     int monte_carlo_steps = 100,
				     string file_name="disks.data")
  {
    int equilibration_steps = monte_carlo_steps / 10;
    cout << " Taking " << equilibration_steps << " equilibration steps" << endl;
    for (int step = 0; step < equilibration_steps; step++)
      monte_carlo_step();
    cout << " Taking " << monte_carlo_steps << " Monte Carlo steps" << endl;
    Matrix<double,1> histogram(number_of_bins);
    for (int step = 0; step < monte_carlo_steps; step++) {
      monte_carlo_step();
      const double pi = 4 * atan(1.0);
      double Delta_A2 = (K * K - 1) * pi * d_0 * d_0 / double(number_of_bins);
      // loop over all pairs of disks
      for (int i = 0; i < N - 1; i++) {
	for (int j = i + 1; j < N; j++) {
	  Matrix<double,1> dr = r[i] - r[j];
	  // if x or y separation > A/2 use closest image separation
	  for (int k = 0; k < 2; k++)
	    if (abs(dr[k]) > 0.5 * A)
	      dr[k] *= 1 - 1 / abs(dr[k]);
	  int bin = int(pi * (dot(dr, dr) - d_0 * d_0) / Delta_A2);
	  if (bin >= 0 && bin < number_of_bins)
	    histogram[bin] += 1;
	}
      }
    }
    ofstream data_file(file_name.c_str());
    for (int bin = 0; bin < number_of_bins; bin++)
      data_file << bin + 1 << '\t' << histogram[bin] << '\n';
    data_file.close();
    cout << " Radial distribution histogram in file " << file_name << endl;
  }

protected : 

  Random rng;                 // random number generator

  double A;                   // side of periodic volume
  int N;                      // number of disk in periodic volume
  Matrix<double,2> r;         // r = (x,y) coordinates of the N disks

  double d_0;                 // disk diameter determined by nu
  double alpha;               // maximum step in x or y for Metropolis algorithm


};

int main(int argc, char *argv[])
{
  cout << " Monte Carlo simulation of hard disk gas\n"
       << " ---------------------------------------\n"
       << " Enter value of nu [0..7]: ";
  double nu;
  cin >> nu;
  Disks disks(nu);
  disks.radial_distribution_histogram(nu);

}
