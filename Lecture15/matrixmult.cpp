#include <iostream>
using namespace std;

#include "linalg.hpp"
using namespace cpt;

int main()
{

  // test the matrix multiplication methods

  Matrix<double,2> A1(2,3), A2(3,2);

  double count = 0.0;
  for ( int i = 0; i < 2; ++i ) {
    for ( int j = 0; j < 3; ++j ) {
      A1[i][j] = count;
      count += 1.0;
    }
  }


  count = 0.0;
  for ( int j = 0; j < 2; ++j ) {
    for ( int i = 0; i < 3; ++i ) {
      A2[i][j] = count;
      count += 1.0;
    }
  }

  Matrix<double,2> A3 = A1*A2;

  cout << "A1 = " << endl << A1 << endl;
  cout << "A2 = " << endl << A2 << endl;
  cout << "A3 = " << endl << A3 << endl;


  Matrix<double,1> v(3);
  v[0] = 10.0;
  v[1] = 0.0;
  v[2] = 0.0;

  
  cout << "v = " << endl << v << endl;
  
  
  Matrix<double,1> prod1 = A1 * v;
  cout << "A1*v = " << endl << prod1 << endl;

  return 0;

}
