#include "cptstd.hpp"
#include "linalg.hpp"
#include "random.hpp"
using namespace cpt;


class Metropolis {

public : 

  Metropolis(int inbins, double iwalker, double idelta, int isteps_to_discard, int isteps_to_skip) : 
    x_walker_0(iwalker), x_walker(iwalker), delta(idelta), steps_to_discard(isteps_to_discard), steps_to_skip(isteps_to_skip),
    bins(400), histogram(bins)
  {
    A[0] = 1.0;
    A[1] = 0.75;
    sigma[0] = 1.0;
    sigma[1] = 0.5;
    center[0] = 0.0;
    center[1] = 6.0;


    // Metropolis
    steps = 0;                   // steps so far
    accepts = 0;                 // steps accepted

    x_max = 10.0;                // maximum extent for plotting
    title = "Metropolis Algorithm for Two Gaussians";
    data_file_name = "metropolis.data";
    gnuplot_file_name = "metropolis.gnu";

  }

  double P (double x)                 // probability distribution
  {
    double p = 0;
    for (int i = 0; i < 2; i++)
      p += A[i]*exp(-(x-center[i])*(x-center[i])/(2*sigma[i]*sigma[i]));
    return p;
  }

  bool metropolis_step()              // return true if step accepted, else false
  {
    double x_trial = x_walker + delta*(2*rng.rand()-1);
    double ratio = P(x_trial) / P(x_walker);

    if (rng.rand() < ratio) {
      x_walker = x_trial;
      return true;
    } else
      return false;
  }


  void take_step()
  {
    // take Monte Carlo steps and bin data
    for (int i = 0; i <= steps_to_skip; i++) {
      monte_carlo_step();
      int bin = int((x_walker+x_max)/(2*x_max)*bins);
      if (bin >= 0 && bin < bins)
	++histogram[bin];
    }

  }

  void write(){
    
    // write data for Gnuplot animation
    ofstream data_file(data_file_name.c_str());
    for (int bin = 0; bin < bins; bin++) {
      double x = -x_max + (bin + 0.5) * (2*x_max)/bins;
      data_file << x << '\t' << histogram[bin] << '\n';
    }
    data_file.close();
  }

  void monte_carlo_step()
  {
    for (int i = 0; i < steps_to_discard; i++)
      metropolis_step();
    if (metropolis_step())
      ++accepts;
    ++steps;
  }

protected : 
  Random rng;

  string title;

  // Gaussians

  double
    A[2],              // amplitudes
    sigma[2],          // widths
    center[2];         // positions of centers


  // Metropolis
  double x_walker_0;             // initial position
  double x_walker;               // current position
  double delta;                  // step size

  int steps_to_discard;           // steps to discard
  unsigned long steps;            // steps so far
  unsigned long accepts;          // steps accepted

  double x_max;                // maximum extent for plotting
  const int bins;               // number of bins in histogram
  Matrix<double,1> histogram;   // histogram bins

  string data_file_name,
    gnuplot_file_name;
  int steps_to_skip;            // when plotting


};

int main () {
  string title;
  double x_walker;
  double delta;
  int monte_carlo_steps;
  int steps_to_discard;
  int steps_to_skip;
  cout << " " << title << endl;
  cout << " Enter starting position x_0: ";
  cin >> x_walker;
  cout << " Enter step size delta: ";
  cin >> delta;
  cout << " Enter number of Monte Carlo steps: ";
  cin >> monte_carlo_steps;
  cout << " Enter number of Metropolis steps per Monte Carlo step: ";
  cin >> steps_to_discard;
  cout << " Enter number of Monte Carlo steps to skip between each plot: ";
  cin >> steps_to_skip;

  int nbins = 400;

  Metropolis metropolis( nbins, x_walker, delta, steps_to_discard, steps_to_skip);

  for ( int i = 0; i < int(monte_carlo_steps / (steps_to_skip + 1)); ++i ) {
    metropolis.take_step();
  }

  metropolis.write();
}
