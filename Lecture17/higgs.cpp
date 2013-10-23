#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "nonlin.hpp"
using namespace cpt;

double V(double x)
{
    return -pow(x, 2.0)/ 2 + pow(x, 4.0)/4;
}

int main()
{
    double acc = 1e-6;
    double guess1 = 0.1, guess2 = -0.1;
    double x_min;
    double V_min = find_minimum< double (*)(double)> (guess1, guess2, V, acc, x_min);
    cout << "V(x) = " << V_min << " at x = " << x_min << endl;

    double x_max;
    double V_max = find_maximum< double (*)(double)> (guess1, guess2, V, acc, x_max);
    cout << "V(x) = " << V_max << " at x = " << x_max << endl;
}
