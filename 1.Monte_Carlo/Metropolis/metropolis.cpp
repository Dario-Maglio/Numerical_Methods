/**
* Implementation of the Metropolis-Hastings algorithm.
*/

//----Preprocessor directives--------------------------------------

#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937 periodo di 2^19937 − 1
#include <cmath>

#define N_STEPS 101000
#define DIM_SAMPLING 1000
#define SEED 42

#define AVERAGE 5.0
#define SIGMA 1.0
#define START 0.
#define DELTA 0.1

using namespace std; 

//----Contents-----------------------------------------------------

/* Probability ratio function */
inline double prf(double q, double q_try) {
    q = pow(q - AVERAGE, 2) - pow(q_try - AVERAGE, 2);
    return exp(q / (2*pow(SIGMA, 2)));
}



/* Main program */
int main() {
    int n=0;
    double x, y, q_try, q=START, sample_average=0., sample_ave_average=0.;
    
    // create files 
    ofstream file_m, file_a;
    file_m.open("f_metropolis.dat");
    file_a.open("f_averages.dat");
    
    // define the PRNG
    mt19937 generator(SEED);
    uniform_real_distribution<double> prng(0.0, 1.0);
    
    // Metropolis-Hastings
    for(int it=0; it!=N_STEPS; it++){
        x = prng(generator);
        y = prng(generator);
        
        q_try = q + DELTA*(1.0 - 2*x);
        x = prf(q, q_try);

        // accept or reject step 
        if (y < x) {                 
            q = q_try;
        }
        
        sample_average = sample_average + q;
        file_m << it << " " << q << endl;
        
        // calculate the average average	
      	if ((it/DIM_SAMPLING) != n) {
	         sample_average = sample_average / DIM_SAMPLING;
	         file_a << n << " " << sample_average << endl;
	         sample_ave_average = sample_ave_average + sample_average;
	         sample_average=0.;
	         n+=1;
         	}
    }
    
    cout << n << " " << sample_ave_average/n << endl; 
    
    // close files
    file_m.close();
    file_a.close();
    
}

