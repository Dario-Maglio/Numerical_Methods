/*******************************************************************************
*
* PRNG using Mersenne and generic Linear congruent generators.
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937
#include <climits>       // import INT_MAX = 2^31 - 1 or 2147483647

#define N 10000          // Mersenne parameters:
#define A 13             // A = 16807       -> a parallel lines
#define C 42             // C = 0
#define M 12345671       // M = INT_MAX     -> affect periodicity
#define SEED 42

using namespace std;

//--- Contents -----------------------------------------------------------------

inline unsigned long int lcg(unsigned long int x) {
    /* Linear congruent generator */
    return (A * x + C) % M;
}

//--- Main ---------------------------------------------------------------------

int main() {
    /* Main program for PRNG */

    unsigned long int x = lcg(SEED);
    double random_lcg = x / (double) M;

    mt19937 generator(SEED);
    uniform_real_distribution<double> dis(0.0, 1.0);
    double random_mars = dis(generator);

    ofstream file;
    file.open ("random.dat");

    for(int it=0; it!=N; it++){
        x = lcg(x);
        random_lcg = (float)x / M;

        random_mars = dis(generator);

        file << it << " " << random_lcg << " " << random_mars << endl;
    }

    file.close();
}
