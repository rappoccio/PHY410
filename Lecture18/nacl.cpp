#include "nacl.hpp"

int main()
{
    ofstream file("ion-potential.data");
    Na na1, na2;
    Cl cl;
    Matrix<double,1> r(3);
    double rmin = 0.01, rmax = 1;
    int points = 100;
    for (int i = 0; i < points; i++) {
        r[0] = rmin + i * (rmax - rmin) / double(points - 1);
        na2.r = r;
        file << r[0] << '\t' << V(na1, na2) << '\n';
    }
    file << '\n';
    for (int i = 0; i < points; i++) {
        r[0] = rmin + i * (rmax - rmin) / double(points - 1);
        cl.r = r;
        file << r[0] << '\t' << V(na1, cl) << '\n';
    }
    file.close();
}
