/**
* PRNG using Mersenne and generic Linear congruent generators.
*/

//----Preprocessor directives--------------------------------------

#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937
#include <climits>       // import INT_MAX = 2^31 - 1 or 2147483647

#define N 10000
#define A 33
#define C 11
#define M INT_MAX
#define SEED 42

using namespace std; 

//----Contents-----------------------------------------------------

/* Linear congruent generator */
inline unsigned int lcg(unsigned int x) {
    return (A * x + C) % M;
}



/* Main program */
int main() {
    unsigned int x = lcg(SEED);
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

