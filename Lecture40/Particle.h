#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include <math.h>
#include <string>
#include <iostream>
#include <sstream>


void getQ_T3_NC(int pdgID, double & q, double & t3, double & NC)
{
  if ( std::abs(pdgID) < 7 ) NC = 3;
  else NC = 1;
  
  if       ( pdgID ==   2 || pdgID ==   4 || pdgID ==   6 ) { q = 2./3.; t3= 0.5;} // up type quark
  else if  ( pdgID ==  -2 || pdgID ==  -4 || pdgID ==  -6 ) { q =-2./3.; t3=-0.5;} // up type anti-quark
  else if  ( pdgID ==   1 || pdgID ==   3 || pdgID ==   5 ) { q =-1./3.; t3=-0.5;} // down type quark
  else if  ( pdgID ==  -1 || pdgID ==  -3 || pdgID ==  -5 ) { q = 1./3.; t3= 0.5;} // down type anti-quark
  else if  ( pdgID ==  11 || pdgID ==  13 || pdgID ==  15 ) { q = -1.0 ; t3=-0.5;} // charged leptons
  else if  ( pdgID == -11 || pdgID == -13 || pdgID == -15 ) { q =  1.0 ; t3= 0.5;} // charged leptons  
  else if  ( pdgID ==  12 || pdgID ==  14 || pdgID ==  16 ) { q =  0.0 ; t3= 0.5;} // neutrinos
  else if  ( pdgID == -12 || pdgID == -14 || pdgID == -16 ) { q =  0.0 ; t3=-0.5;} // anti-neutrinos  
  else {
    std::cout << "Error with PDG ID " << pdgID << std::endl;
    q = -9999; t3= -9999;
  }  
}

class Particle{
public :
  Particle( LorentzVector const & p4, int pdgID ) :
    p4_(p4), pdgID_(pdgID)
  {
    double NC=0;
    getQ_T3_NC(pdgID, q_, t3_, NC);
  }

  inline const  LorentzVector & p4() const { return p4_;}
  inline double q() const { return q_;}
  inline double t3() const { return t3_;}
  inline int pdgID() const { return pdgID_; }
  
private :
  LorentzVector p4_; // 4-vector
  int pdgID_;     // PDG ID from http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
  double q_;      // electric charge
  double t3_;     // weak isospin

};

#endif
