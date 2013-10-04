#include <iostream>
using namespace std;

#include "linalg.hpp"
using namespace cpt;

int main()
{


    double m1 = 1, m2 = 2, m3= 1;

    // Matrix with masses as the diagonal elements
    Matrix<double,2> M(3, 3);
    M(0,0) = m1;
    M(1,1) = m2;
    M(2,2) = m3;
    cout << "M =\n" << M;


    Matrix<double,2> Minv(3,3);
    inverse(M, Minv);
    cout << "Minv=" << endl << Minv << endl;

    // "Spring" constants affecting each mass
    double k12 = 1, k23 = 1;
    double Lagrange[3][3] = {
        {   k12,    -k12,       0       },
        {   -k12,   k12 + k23,  -k23    },
        {   0,      -k23,       k23     }
    };
    Matrix<double,2> K(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            K(i,j) = Lagrange[i][j];
    cout << "K =\n" << K;




    // Solve with generalized eigenvector solution
    Matrix<double,1> eigenvalues = solve_eigen_generalized(K, M, true);

    cout << "Eigenvalues =\n" << eigenvalues << endl
         << "Eigenvectors =\n" << K;


}
