// Adapted from the MSTW2008 code file example.cc
// Needs mstwpdf.cc mstwpdf.h mstw2008lo.00.dat
// c++  mstw.cpp  mstwpdf.cc

#include <cmath>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>

#include "mstwpdf.h"
#include "LorentzVector.h"
#include "Particle.h"


// Calculate Drell-Yan cross section given initial quark momenta p1 and p2
double dSigmaDOmega( const Particle & p1, const Particle & p2, int pdgID_out, double cosTheta ) {
  if ( p1.pdgID() != -p2.pdgID() ) {
    std::cout << "Flavor changing neutral currents are not observed in the SM." << std::endl;
    return 0;
  }
  // Constants from http://pdg.lbl.gov/2014/reviews/rpp2014-rev-phys-constants.pdf
  static const double MZ = 91.1876;      // GeV
  static const double GF = 1.1663787e-5; // 1/GeV^2
  static const double sin2ThetaW = 0.23126;
  static const double Gamma_Z = 2.4952;

  double t3 = p1.t3();              // z-component of weak isospin of incoming particles
  double q = p1.q() ;               // charge of incoming particles
  auto tot = p1.p4() + p2.p4();     // p4 for the COM collision
  double s = tot.mass2();           // Mandelstam s (sqrt(s) = COM energy)
  double Nc_o = 1.;                 // Number of colors of outgoing particles
  double Nc_i = 1.;                 // Number of colors of incoming particles
  double q_o = -1.;                 // Charge of outgoing particles
  double t3_o = 0.;                 // Weak isospin of outgoing particles  
  getQ_T3_NC(pdgID_out, q_o, t3_o, Nc_o);

  // Compute some constants
  double c2 = 8*GF*MZ*MZ / sqrt(2);
  double L =    ( t3-sin2ThetaW*q);
  double R = -1*(sin2ThetaW*q);
  double Lp =  ( t3_o-sin2ThetaW*q_o);
  double Rp = -1*(sin2ThetaW*q_o);
  // L, R, Lp, Rp are correct
  double dsigma_domega =  Nc_o / ( Nc_i ) * s /( 256. * M_PI*M_PI) / ( pow( s-MZ*MZ, 2.0) + s* Gamma_Z*Gamma_Z) *
    c2*c2* ( (L*L+R*R)*(Lp*Lp+Rp*Rp)*(1+cosTheta*cosTheta) + (L*L-R*R)*(Lp*Lp-Rp*Rp)*2*cosTheta ) ;

  return dsigma_domega;
}



int main()
{

  // create incoming e+ and e-
  Particle em_i   ( LorentzVector(0,0, 91.2*0.5,91.2*0.5),  11);
  Particle ep_i   ( LorentzVector(0,0,-91.2*0.5,91.2*0.5), -11);

  // step through cos(theta) values from -1 to 1
  // create an output plot for use in matplotlib
  std::ofstream outplot("ee_Z_ee_dsigma_domega.txt");
  double sigma = 0.;
  double dcosTheta = 0.05;
  for ( double cosTheta = -1.0; cosTheta < 1.0; cosTheta += dcosTheta ) {
    // Produce e+e- --> Z --> e+e-
    double ds_dO = dSigmaDOmega( em_i, ep_i, 11, cosTheta );
    char buf[1000];
    sprintf(buf, "%6.2f %10.6e", cosTheta, ds_dO);
    outplot << buf << std::endl;
    // Integral = sum( dsigma/domega) * dcos(theta)
    sigma += ds_dO * dcosTheta;
  }
  std::cout << "Sigma = " << sigma << std::endl;


  // Print the PDFs
  // Get the PDF file
  string grid_file_name("Grids/mstw2008lo.00.dat");
  c_mstwpdf *pdf = new c_mstwpdf(grid_file_name.c_str());


  // Here the PDG notation is used for the parton flavour
  // (apart from the gluon has f=0, not 21):
  //  f =   -6,  -5,  -4,  -3,  -2,  -1,0,1,2,3,4,5,6
  //    = tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t.
  // Can also get valence quarks directly:
  //  f =  7, 8, 9,10,11,12.
  //    = dv,uv,sv,cv,bv,tv.
  // Photon: f = 13.
  // 

  std::ofstream pdfout("proton_pdfs.txt");
  // Specify the momentum fraction "x" and scale "q".
  double q = 100.;

  // Print out some values
  cout << "q = " << q << ", "
       << "alphaS(Q0) = " << pdf->alphaSQ0 << ", "
       << "alphaS(MZ) = " << pdf->alphaSMZ << ", "
       << "alphaSorder = " << pdf->alphaSorder << ", "
       << "alphaSnfmax = " << pdf->alphaSnfmax << endl;

  // Just to get a reasonable sampling, do logarithmic binning, then progress to linear binning
  double xvals[] = { 0.001, 0.002, 0.004, 0.008, 0.01, 0.02, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
  for ( auto x : xvals ) {    
    pdf->update(x,q);
    double upv = x * pdf->cont.upv;
    double dnv = x * pdf->cont.dnv;
    double usea = x * pdf->cont.usea;
    double dsea = x * pdf->cont.dsea;
    double str = x * pdf->cont.str;
    double sbar = x * pdf->cont.sbar;
    double chm = x * pdf->cont.chm;
    double cbar = x * pdf->cont.cbar;
    double bot = x * pdf->cont.bot;
    double bbar = x * pdf->cont.bbar;
    double glu = x * pdf->cont.glu;
    double phot = x * pdf->cont.phot;

    char buf[1000];
    sprintf( buf, "%8.4e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e", x, upv, dnv, usea, dsea, str, sbar, glu );
    pdfout << buf << std::endl;
  }

  // Print out grid ranges, heavy quark masses, and alphaS values


  delete pdf;
  pdfout.close();
  outplot.close();
  return 0;
}
