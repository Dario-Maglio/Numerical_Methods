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
#define N_STEPS 10015000           // dim of sample to generate
#define DIM_SAMPLING 10000         // dim of subsamples for calculating averages
#define DIM_SUBSAMP 4194304        // dim of subsample for autocorrelation 2^22
#define DIM_THERM 15000            // dim of subsample for thermalization
#define NUM_FAKESAMPLES 100        // num of fake samples

#define AVERAGE 5.0
#define SIGMA 1.
#define START 0.
#define DELTA 0.1                  // delta circa 3-4 sigma gives the best err

using namespace std;

// define the PRNG
mt19937_64 generator(SEED);
uniform_real_distribution<double> prng(0.0, 1.0);
uniform_int_distribution randomint(0, DIM_SUBSAMP);

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
double load_and_average(int dim_therm, int dim_subsample, vector<double>& x){
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
double autocorrelation(vector<double>& x, float x_ave){
    int k = 0, steps = size(x);
    double value = 1., tau_int = 0., var_bias = 0., sigma;

    // compute the biased variance
    for (int j = 0; j < steps; j++) var_bias += pow(x[j] - x_ave , 2);
    var_bias = var_bias / (steps - 1);
    cout << "logging: biased std dev is " << sqrt(var_bias/steps) << endl;

    // compute tau_int and the autocorrelation function
    ofstream file;
    file.open("file_autocor.dat");

    while (value > 0.) {
        value = 0.;
        for (int i = 0; i < (steps - k); i++){
            value += pow(x_ave, 2) + (x[i] * x[i+k]) - x_ave * (x[i] + x[i+k]);
        }
        value = value / (steps - k);
        tau_int += value;
        k++;

        file << k << " " << value << " " << tau_int << endl;
        cout << "logging: C at " << k << " is " << value << endl;
    }
    file.close();

    sigma = (1 + 2*tau_int) / steps;
    sigma = sqrt(var_bias * sigma);
    return sigma;
}

/* Bootstrap without autocorrelation */
double bootstrap_step(vector<double>& x){
    int steps = size(x);
    double x_fake, fake_ave = 0.;
    vector<double> fake;
    fake.reserve(steps);

    for (int i = 0; i < steps; i++){
        x_fake = x[randomint(generator)];
        fake_ave += x_fake;
        fake.push_back(x_fake);
    }
    fake.clear();

    fake_ave = fake_ave / steps;
    return fake_ave;
}

double bootstrap(vector<double>& x, int num_fake_samples){
    double fake_ave, fake_ave_ave = 0., sigma = 0.;
    vector<double> fake_averages;
    fake_averages.reserve(num_fake_samples);

    ofstream file_boot;
    file_boot.open("file_boot.dat");

    for (int k = 5; k <= num_fake_samples; k += 5){
        fake_ave_ave = 0.;
        for (int i = 0; i < k; i++){
            fake_ave = bootstrap_step(x);
            fake_ave_ave += fake_ave;
            fake_averages.push_back(fake_ave);
        }
        fake_ave_ave = fake_ave_ave / k;

        sigma = 0.;
        for (int i = 0; i < k; i++){
            sigma += pow(fake_averages[i], 2);
        }
        sigma = sigma / k;
        sigma = sigma - pow(fake_ave_ave, 2);
        sigma = sqrt(sigma);

        fake_averages.clear();
        file_boot << k << " " << sigma << endl;
        cout << "logging: sigma at " << k << " is " << sigma << endl;
    }

    file_boot.close();
    return sigma;
}

/* Bootstrap with autocorrelation */
double bootstrap_corr_step(vector<double>& x, int k){
    int index, num_block, steps = size(x);
    double x_fake, fake_ave = 0.;
    vector<double> fake;
    fake.reserve(steps);

    num_block = steps / pow(2, k);

    for (int i = 0; i < num_block; i++){
        index = randomint(generator);
        for (int j = 0; j < pow(2, k); j++){
            if (index + j == steps) index = index - steps;
            x_fake = x[index + j];
            fake_ave += x_fake;
            fake.push_back(x_fake);
        }
    }
    fake.clear();

    fake_ave = fake_ave / steps;
    return fake_ave;
}

double bootstrap_corr(vector<double>& x, int num_fake_samples){
    double fake_ave, fake_ave_ave = 0., sigma = 0.;
    vector<double> fake_averages;
    fake_averages.reserve(num_fake_samples);

    ofstream file_boot;
    file_boot.open("file_boot_corr.dat");

    for (int k = 5; k < 16; k++){
        fake_ave_ave = 0.;
        for (int i = 0; i < num_fake_samples; i++){
            fake_ave = bootstrap_corr_step(x, k);
            fake_ave_ave += fake_ave;
            fake_averages.push_back(fake_ave);
        }
        fake_ave_ave = fake_ave_ave / num_fake_samples;

        sigma = 0.;
        for (int i = 0; i < num_fake_samples; i++){
            sigma += pow(fake_averages[i], 2);
        }
        sigma = sigma / num_fake_samples;
        sigma = sigma - pow(fake_ave_ave, 2);
        sigma = sqrt(sigma);

        fake_averages.clear();
        file_boot << k << " " << sigma << endl;
        cout << "logging: sigma at " << k << " is " << sigma << endl;
    }

    file_boot.close();
    return sigma;
}

//---- Main --------------------------------------------------------------------

int main() {
    double x_ave, sigma;
    vector<double> data;

    // generate_MC_data();

    x_ave = load_and_average(DIM_THERM, DIM_SUBSAMP, data);
    cout << "The average of the sample is " << x_ave << endl;

    // sigma = autocorrelation(data, x_ave);
    // cout << "The result is " << x_ave << " +- " << sigma << endl;

    sigma = bootstrap(data, NUM_FAKESAMPLES);
    cout << "Bootstrap gives: " << sigma << endl;

    sigma = bootstrap_corr(data, NUM_FAKESAMPLES);
    cout << "Bootstrap corr gives: " << sigma << endl;

    return 0;
}
