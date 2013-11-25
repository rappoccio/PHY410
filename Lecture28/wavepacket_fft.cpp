#include "wavepacket.hpp"

class WavepacketFFT : public WavepacketBase {
public : 

 WavepacketFFT( int iN=128, double iL=100, double idt=0.1, bool iper=true ) : 
   WavepacketBase( iN, iL, idt, iper ),
    T_exp_factor(iN),
    V_exp_factor(iN)
  {
    for (int j = 0; j < N; j++) {
      // kinetic factor exp[-iT/h_bar dt]
      const double pi = 4 * atan(1.0);
      double p = j < N / 2 ? j : j - N;
      p *= h_bar * 2 * pi / L;
      double theta = - p * p / (2 * mass) / h_bar * dt;
      T_exp_factor[j] = complex<double>(cos(theta), sin(theta));

      // potential factor exp[-iV(x)/(2h_bar) dt]
      theta = - V(x[j]) / 2 / h_bar * dt;
      V_exp_factor[j] = complex<double>(cos(theta), sin(theta));
    }
  }

  void take_step() {

    // first half potential phase rotation
    for (int j = 0; j < N; j++)
      psi[j] *= V_exp_factor[j];

    // FFT to momentum space
    fft(psi);

    // kinetic phase rotation
    for (int j = 0; j < N; j++)
      psi[j] *= T_exp_factor[j];

    // FFT back to position space
    bool do_inverse = true;
    fft_inv(psi);

    // second half potential phase rotation
    for (int j = 0; j < N; j++)
      psi[j] *= V_exp_factor[j];

    t += dt;
  }

protected : 
  Matrix<complex<double>,1>
  T_exp_factor, V_exp_factor;   // precomputed phase rotations


};



int main()
{


  cout << " Quantum Wavepacket Motion via FFT" << endl;
  WavepacketFFT wavepacket;
  cout << "Got a wavepacket" << endl;

  ofstream file("potential.data");
  for (int i = 0; i < wavepacket.get_N(); i++)
    file << wavepacket.get_x(i) << '\t' << wavepacket.V(wavepacket.get_x(i) ) << '\n';
  file.close();
  cout << " saved V(x) in file potential.data" << endl;



#ifdef _WIN32
  ostringstream pyout;
  pyout << "env python.exe animator_for_cpp.py " << wavepacket.get_N();
  FILE *pypipe = _popen(pyout.str().c_str(), "w");
#else
  ostringstream pyout;
  pyout << "/usr/bin/env python animator_for_cpp.py " << wavepacket.get_N();
  FILE *pypipe = popen(pyout.str().c_str(), "w");
#endif


  wavepacket.save_psi(0);
  int plots = 10;
  for (int plot = 1; plot <= plots; plot++) {
    cout << "plot = " << plot << endl;
    double delta_t = 0;
    while (delta_t < wavepacket.get_L() / (plots * wavepacket.get_velocity())) {
      wavepacket.take_step();
      delta_t += wavepacket.get_dt();
    }
    wavepacket.save_psi(plot);
  }

  
  // simple Mpl pipes animation
  cout << " Enter animation time: ";
  double t_max;
  cin >> t_max;
  wavepacket = WavepacketFFT();
  double frame_rate = 30;
  double dt_frame = 1 / frame_rate;
  int steps_per_frame = max(1, int(dt_frame /   wavepacket.get_dt()));

  while (wavepacket.get_t() < t_max) {
    ostringstream os;
    for (int i = 0; i < wavepacket.get_N(); i++)
      os << std::abs(wavepacket.get_psi()[i]) << ',';
    fprintf(pypipe, "%s\n", os.str().c_str());
    fflush(pypipe);
    time_t start_time = clock();
    for (int step = 0; step < steps_per_frame; step++) {
      wavepacket.take_step();
    }
    while (true) {
      double secs = (clock() - start_time) / double(CLOCKS_PER_SEC);
      if (secs > dt_frame)
	break;
    }
  }
  fclose(pypipe);
}

