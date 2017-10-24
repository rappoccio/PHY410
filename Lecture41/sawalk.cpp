#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

inline double std_rand()
{
    return rand() / (RAND_MAX + 1.0);
}

struct Site {               // object to represent lattice site
    int x;                  // x coordinate
    int y;                  // y coordinate
};

int main() {

    cout << " Self-Avoiding Walks on a Square Lattice\n"
         << " ---------------------------------------\n"
         << " Enter number of steps in walk: ";
    int n_steps, n_walks;
    cin >> n_steps;
    cout << " Enter number of walks to generate: ";
    cin >> n_walks;

    int walks = 0;
    int failed_walks = 0;
    double r2av = 0;
    double r4av = 0;

    // generate walks
    while (walks < n_walks) {
        vector<Site> sites;         // set of occupied lattice sites
        Site s;
        s.x = 0;
        s.y = 0;
        sites.push_back(s);
        bool walk_failed = false;

        // loop over desired number of steps
        for (int step = 0; step < n_steps; step++) {

            // take a random step
            double d = std_rand();
            if      (d < 0.25)  ++s.x;   // step East
            else if (d < 0.50)  ++s.y;   // step North
            else if (d < 0.75)  --s.x;   // step West
            else                --s.y;   // step South

            // check whether the site is occupied
            bool occupied = false;
            for (int i = 0; i < sites.size(); i++) {
                if (s.x == sites[i].x && s.y == sites[i].y) {
                    occupied = true;
                    break;
                }
            }
            if (occupied) {
                walk_failed = true;
                break;
            }

            sites.push_back(s);
        }

        if (walk_failed) {
            ++failed_walks;
            continue;
        }

        double r2 = s.x * s.x + s.y * s.y;
        r2av += r2;
        r4av += r2 * r2;
        ++walks;
    }

    r2av /= n_walks;
    r4av /= n_walks;
    double stdDev = sqrt(r4av - r2av * r2av);
    double totalWalks = n_walks + failed_walks;
    double failedPercent = failed_walks / totalWalks * 100.0;
    cout << " Mean square distance <r^2> = " << r2av << "\n"
         << " Standard deviation         = " << stdDev << "\n"
         << " Percentage failed walks    = " << failedPercent
         << endl;
}
