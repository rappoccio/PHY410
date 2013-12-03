#include "cptstd.hpp"
#include "random.hpp"
using namespace cpt;

int poor_seed = 123;                    // maintain seed for poor generator
double poor() {                         // example from numerical recipes
    const int IM = 6075, IA = 106, IC = 1283;
    int seed = poor_seed;
    seed = (IA * seed + IC) % IM;       // linear congruential algorithm
    poor_seed = seed;
    return seed / double(IM);           // make range 0.0 ... 1.0
}

int park_miller_seed = 1;
double park_miller() {                  // Park-Miller generator
    const int IA = 16807, IC = 2147483647, IQ = 127773, IR = 2836;
    int seed = park_miller_seed;
    int h = seed / IQ;
    int l = seed % IQ;
    seed = IA * l - IR * h;
    if (seed <= 0)
        seed += IC;
    park_miller_seed = seed;
    return seed / double(IC);
}

int std_rand_seed;
double std_rand() {                     // C++ standard library generator
    std_rand_seed = rand();
    return std_rand_seed / (RAND_MAX + 1.0);
}

Random rg;
double xorshift() {
    double r = rg.rand();
    return r;
}

double mersenne() {
    double r = rg.rand();
    return r;
}

void test(double generator(), string name, string remark) {
    cout << " " << name << "() -- " << remark << endl;

    cout << " checking period ... " << flush;
    double start = generator();
    long int i, max_steps = 10000000L;
    for (i = 0L; i < max_steps; i++) {
        if (generator() == start) {
            cout << " repeats after "  << i << " steps\n";
            break;
        }
    }
    if (i >= max_steps)
        cout << " period larger than " << max_steps << endl;

    clock_t start_time = clock();
    const int BINS = 10000;
    cout << " binning " << max_steps << " tries in "
         << BINS << " bins ..." << flush;
    int bin[BINS];
    for (int b = 0; b < BINS; b++)
        bin[b] = 0;
    for (i = 0L; i < max_steps; i++) {
        int b = int(generator() * BINS);
        if (b == BINS)
          --b;
        ++bin[b];
    }
    double chiSqr = 0;
    double expect = max_steps / double(BINS);
    for (int b = 0; b < BINS; b++) {
        double diff = bin[b] - expect;
        chiSqr += diff * diff;
    }
    chiSqr /= expect;
    cout << " chi-square/d.o.f. = " << chiSqr/(BINS-1) << endl;
    cout << " CPU time = " << (clock() - start_time) / double(CLOCKS_PER_SEC)
         << " seconds\n";

    string fileName = name + ".data";
    cout << " writing eyeball test data to " << fileName << endl;
    ofstream dataFile(fileName.c_str());
    for (int i = 0; i < 10000; i++)
        dataFile << generator() << '\t' << generator() << '\n';
    dataFile.close();
    cout << endl;
}

int main() {
    cout << " Testing random number generators\n"
         << " --------------------------------\n\n";
    test(poor, "poor", "from Numerical Recipes period 6,075");
    test(park_miller, "park_miller", "Park-Miller multiplicative generator");
    test(std_rand, "std_rand", "rand() from C++ Standard Library");
    rg.set_xorshift();
    test(xorshift, "xorshift", "Marsaglia Xorshift generator from random.hpp");
    rg.set_mersenne_twister();
    test(mersenne, "mersenne", "Mersenne Twister generator");
}
