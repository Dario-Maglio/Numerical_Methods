/*******************************************************************************
*
* Implementation of the Metropolis-Hastings algorithm.
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------

#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937 periodo di 2^19937 − 1
#include <cmath>
#include <vector>

#define N_STEPS 1015000
#define DIM_SAMPLING 1000
#define DIM_SUBSAMP 1000000
#define DIM_THERM 10000
#define SEED 42

#define AVERAGE 5.0
#define SIGMA 1.
#define START 0.
#define DELTA 0.1

using namespace std;

//---- Contents ----------------------------------------------------------------

/* Probability ratio function */
inline double prf(double q, double q_try) {
    q = pow(q - AVERAGE, 2) - pow(q_try - AVERAGE, 2);
    return exp(q / (2*pow(SIGMA, 2)));
}

/* Autocorrelation function, variance and error */
void auto_corre(int dim_subsample){
    int steps = 0;
    double value, x_ave = 0., tau_int = 0., var_bias = 0., sigma;
    vector<double> x;

    // load the data file and compute the average
    ifstream data("f_metropolis.dat");
    if (data.is_open()) {
        while ((data >> value) && (steps < DIM_THERM)) steps += 1;
        steps = 0;
        while ((data >> value) && (steps <= dim_subsample)){
            steps += 1;
            x_ave += value;
            x.push_back(value);
        }
        data.close();
        x_ave = x_ave / steps;
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
    cout << "The average of the sample is " << x_ave << endl;

    // compute the biased variance
    for (int k = 0; k < steps; k++) var_bias += pow(x[k] - x_ave , 2);
    var_bias = (var_bias/steps)  / (steps - 1);
    cout << "The biased std deviation is " << sqrt(var_bias) << endl;

    // compute tau_int and the autocorrelation function
    ofstream file_corr, file_tau;
    file_corr.open("f_autocorr.dat");
    file_tau.open("f_tau_int.dat");

    for (int k = 0; k < steps; k++){
        value = 0.;
        for (int i = 0; i < (steps - k); i++){
            value += pow(x_ave, 2) + (x[i] * x[i+k]) - x_ave * (x[i] + x[i+k]);
        }
        value = value / (steps - k);
        file_corr << value << endl;
        tau_int += value;
        file_tau << tau_int << endl;

        // the approximation starts breaking down; steps must be >> t_int
        if (value < 0) break;
    }

    file_corr.close();
    file_tau.close();

    sigma = sqrt(var_bias * (1 + 2*tau_int));
    cout << "The result is " << x_ave << " +- " << sigma << endl;
}

/* generate the MC data */
void generate_MC_data() {
    double x, y, q_try, q=START, sample_average=0.;

    // create files
    ofstream file_m, file_a;
    file_m.open("f_metropolis.dat");
    file_a.open("f_averages.dat");

    // define the PRNG
    mt19937 generator(SEED);
    uniform_real_distribution<double> prng(0.0, 1.0);

    // Metropolis-Hastings
    for(int it=1; it<=N_STEPS; it++){
        x = prng(generator);
        y = prng(generator);

        q_try = q + DELTA*(1.0 - 2*x);
        x = prf(q, q_try);

        // accept or reject step
        if (y < x) {
            q = q_try;
        }

        sample_average = sample_average + q;
        file_m << q << endl;

        // compute the average over a sample
      	if ((it % DIM_SAMPLING) == 0) {
	         sample_average = sample_average / DIM_SAMPLING;
	         file_a << sample_average << endl;
	         sample_average=0.;
         	}
    }

    // close files
    file_m.close();
    file_a.close();
}

//---- Main --------------------------------------------------------------------

int main() {
    // generate_MC_data();

    auto_corre(DIM_SUBSAMP);

    return 0;
}
