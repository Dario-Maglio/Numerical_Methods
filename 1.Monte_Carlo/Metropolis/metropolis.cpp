/*******************************************************************************
*
* Implementation of the Metropolis-Hastings algorithm with error estimation.
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>

#define SEED 42
#define N_STEPS 1015000            // dim of sample to generate
#define DIM_SAMPLING 1000          // dim of subsamples for calculating averages
#define DIM_SUBSAMP 1000000        // dim of subsample for autocorrelation
#define DIM_THERM 10000            // dim of subsample for thermalization

#define AVERAGE 5.0
#define SIGMA 1.
#define START 0.
#define DELTA 0.1                  // delta circa 3-4 sigma gives the best err

using namespace std;

//---- Contents ----------------------------------------------------------------

/* Probability ratio function */
inline double prf(double q, double q_try) {
    q = pow(q - AVERAGE, 2) - pow(q_try - AVERAGE, 2);
    return exp(q / (2*pow(SIGMA, 2)));
}

/* generate the MC data */
void generate_MC_data() {
    double x, y, q_try, q=START, sample_average=0.;

    // create files
    ofstream file_m, file_a;
    file_m.open("file_metropol.dat");
    file_a.open("file_averages.dat");

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

/* Load the sample data and compute the average */
float load_and_average(int dim_therm, int dim_subsample, vector<double>& x){
  int steps = 0;
  double value, x_ave = 0.;

  // load the data file and compute the average
  ifstream data("file_metropol.dat");
  if (data.is_open()) {
      while ((data >> value) && (steps < dim_therm)) steps += 1;
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

  return x_ave;
}

/* Autocorrelation function, variance and error */
float auto_corre(vector<double>& x, float x_ave){
    int steps = size(x);
    double value, tau_int = 0., var_bias = 0., sigma;

    // compute the biased variance
    for (int k = 0; k < steps; k++) var_bias += pow(x[k] - x_ave , 2);
    var_bias = var_bias / (steps - 1);
    cout << "logging: biased std dev is " << sqrt(var_bias/steps) << endl;

    // compute tau_int and the autocorrelation function
    ofstream file_cor, file_tau;
    file_cor.open("file_autocor.dat");
    file_tau.open("file_tau_int.dat");

    for (int k = 0; k < steps; k++){
        value = 0.;
        for (int i = 0; i < (steps - k); i++){
            value += pow(x_ave, 2) + (x[i] * x[i+k]) - x_ave * (x[i] + x[i+k]);
        }
        value = value / (steps - k);
        file_cor << value << endl;
        tau_int += value;
        file_tau << tau_int << endl;

        // the approximation starts breaking down; steps must be >> t_int
        if (value < 0) break;
    }

    file_cor.close();
    file_tau.close();

    sigma = (1 + 2*tau_int) / steps;
    sigma = sqrt(var_bias * sigma);
    return sigma;
}

//---- Main --------------------------------------------------------------------

int main() {
    float x_ave, sigma;
    vector<double> data;

    generate_MC_data();

    x_ave = load_and_average(DIM_THERM, DIM_SUBSAMP, data);
    cout << "The average of the sample is " << x_ave << endl;

    sigma = auto_corre(data, x_ave);
    cout << "The result is " << x_ave << " +- " << sigma << endl;

    return 0;
}
