#include <cmath>
#include <iostream>
using namespace std;
 
#include "linalg.hpp"
using namespace cpt;
 
int main()
{
    cout << " Unbalanced Wheatstone bridge equations\n"
         << " --------------------------------------\n";
 
    double v0 = 1.5,
           r1 = 100, r2 = r1, r3 = 150,
           rx = 120, ra = 1000, rv = 10;
 
    Matrix<double,2> v(3, 1);       // column vector with 3 rows
    v[0][0] = v0;
    cout << " v = \n" << v << endl;
    cout << "v.dim1 = " << v.dim1() << endl;
 
    Matrix<double,2> R(3, 3);       // 3x3 resistance matrix
    R[0][0] = r1 + rv;              // set components using slicing notation
    R[0][1] = r2;
    R[0][2] = rv;
 
    R[1][0] = r1 + ra;
    R[1][1] = -ra;
    R[1][2] = -r3;
 
    R(2,0) = rx + ra;               // or use subscripting notation
    R(2,1) = -r2 - rx - ra;
    R(2,2) = rx;
 
    cout << " R = \n" << R << std::endl;
    cout << " R.dim1 = " << R.dim1() << endl;
    cout << " R.dim2 = " << R.dim2() << endl;
 
    // the solve_Gauss_Jordan replaces R by R^-1 and v by i
    // so save the original R and copy v into a vector i
    Matrix<double,2> R_save(R), i(v);
 
    solve_Gauss_Jordan(R, i);
    cout << " Solution using Gauss-Jordan elimination\n i = \n"
         << i << endl;
 
    // find the other currents
    cout << " i_a = i_1 - i_2 = " << i[0] - i[1] << '\n'
         << " i_v = i_1 + i_3 = " << i[0] + i[2] << '\n'
         << " i_x = i_1 - i_2 + i_3 = " << i[0] - i[1] + i[2] << '\n'
         << endl;
 
    // see whether LU decomposition gives the same answer
    i = v;
    solve_LU_decompose(R_save, i);
    cout << " Solution using LU Decompositon\n i = \n"
         << i << endl;
}
