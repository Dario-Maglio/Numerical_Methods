/**
* Implementation of the Metropolis algorithm.
*/

//----Preprocessor directives--------------------------------------

#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937
#include <math.h>

#define N 10000
#define SEED 42

#define AVERAGE 5.0
#define SIGMA 1.0
#define START 0.
#define DELTA 0.1

using namespace std; 

//----Contents-----------------------------------------------------

/* Probability ratio function */
inline double pdf(double q, double q_try) {
    return exp(((q - AVERAGE)**2 - (q_try - AVERAGE)**2)/2.0/(SIGMA*SIGMA));
}



/* Main program */
int main() {
    int j=0, k=0, n=0;
    double x, y, q_try, q=START, sample_average=0., sample_ave_average=0.;
    
    // create files 
    ofstream file_m, file_a;
    file_m.open("f_metropolis.dat");
    file_a.open("f_averages.dat");
    
    // define the PRNG
    mt19937 generator(SEED);
    uniform_real_distribution<double> prng(0.0, 1.0);
    
    // initial configuration
    
    for(int it=0; it!=N; it++){
        x = pnrg(generator);
        y = pnrg(generator);
        
        q_try = q + delta*(1.0 - 2.0*x);
    
        x = fun(q, q_try);

        if (y.lt.z) then                    //accept reject 
            q = q_try
            acc = 1.0                       //accettanza = 1 o 0
 	    j=j+1
 	    k=k+1
 	    sampave = sampave + q
            write(1,*) j,q
        else
            acc = 0.0                       //se non accetto tengo q vecchio
        endif                 
	 
	if ((i/10000) /= n) then
	    sampave = sampave / k
	    write(2,*) n, sampave
	    sampaveave = sampaveave + sampave
	    n = n + 1
	    sampave=0.
	    k=0
	endif
        file << it << " " << random_lcg << " " << random_mars << endl;
    }
    
    write(*,*) sampaveave/n
    
    // close files
    file_m.close();
    file_a.close();
    
}

