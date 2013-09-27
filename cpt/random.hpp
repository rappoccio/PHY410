#ifndef CPT_RANDOM_HPP
#define CPT_RANDOM_HPP

#include <stdint.h>
// define a common name Ullong for unsigned 64-bit integers
#ifdef _MSC_VER   /* Microsoft C++ */
  typedef unsigned __int64 Ullong;    // 64-bit unsigned integer
#else             /* Macintosh, Linux */
  typedef unsigned long long Ullong;
#endif

#include <string>

namespace cpt
{

    // This class implements the Park-Miller algorithm,
    // 32- and 64-bit Mersenne Twister algorithms by Nishimura and Matsumoto,
    // & 2 Marsaglia Xor-shift algorithms given in Numerical Recipes 3rd Edition

    class Random {
    public:

        Random();                       // default constructor

        double rand();                  // uniform deviate in [0,1]

        double rand_gauss(              // Gaussian deviate
            double center=0.0,          // default centered at x = 0
            double width=1.0);          // default width sigma = 1.0

        std::string get_algorithm();    // algorithm currently in use

        uint64_t get_seed();            // current value of seed

        void set_park_miller();         // Park-Miller algorithm

        void set_mersenne_twister();    // 32-bit Mersenne twister algorithm

        void set_mersenne_twister_64(); // 64-bit Mersenne twister algorithm

        void set_xorshift();            // Marsaglia Xor-shift period ~10^19

        void set_xorshift_long();       // Marsaglia Xor-shift period ~10^57

        void set_seed(                  // set current seed
            const uint64_t);

        void set_seed_time();           // set seed using std::clock()

    private:

        int algo;
        uint32_t seed_32;
        uint64_t seed_64, u, v, w;

        void init();
        void step();

    };

}   /* namespace cpt */

#endif /* CPT_RANDOM_HPP */
