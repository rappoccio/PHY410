#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

double f(double x) {
    return exp(x) * log(x) - x * x;
}

void print(int step, double x, double dx) {
    cout.setf(ios::right, ios::adjustfield);
    cout << " " << setw(4) << step << "    ";
    cout.precision(15);
    cout.setf(ios::left, ios::adjustfield);
    cout.setf(ios::showpoint | ios::fixed);
    cout << setw(20) << x << "  " << setw(20) << dx << '\n';
}

int main() {

    cout << " Simple search for root of exp(x)*log(x) - x*x\n"
         << " ---------------------------------------------\n"
         << " Enter guess x, step dx, and desired accuracy: ";
    double x, dx, acc;
    cin >> x >> dx >> acc;

    cout << " Step            x                    dx\n"
         << " ----    ------------------    ------------------\n";
    int step = 0;
    print(step, x, dx);
    double f_old = f(x);
    while (fabs(dx) > abs(acc)) {
        x += dx;
        if ( f_old * f(x) < 0 ) {
            x -= dx;
            dx /= 2;
        }
        ++step;
        print(step, x, dx);
    }
}
