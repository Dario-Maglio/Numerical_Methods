/*******************************************************************************
*
* Implementation of the Metropolis-Hastings algorithm with error estimation.
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>

#define SEED 42
#define N_STEPS 10015000           // dim of sample to generate
#define DIM_SUBSAMP 1000000        // dim of subsample to study 2^22 4194304
#define DIM_THERM 15000            // dim of subsample for thermalization

#define DIM_SAMPLING 10000         // dim of subsamples for calculating averages
#define NUM_FAKESAMPLES 300        // num of fake samples bootstrap
#define CORREL_LENGHT 65536        // dim of the correlated block bootstrap 2^16

#define AVERAGE 0.                 // 5. in other cases
#define SIGMA 1.
#define START 0.
#define DELTA 3.                   // 0.1 in other cases

using namespace std;

constexpr double pi = 3.14159265358979323846;

// Define the PRNG
mt19937_64 generator(SEED);
uniform_real_distribution<double> prng(0.0, 1.0);
uniform_int_distribution<long int> randomint(0, DIM_SUBSAMP);

//---- Contents ----------------------------------------------------------------

void generate_BM_data(){
    /* Generate random numbers with the Box-Muller algorithm */

    double x1, x2, y1, y2;
    ofstream file_b;
    file_b.open("file_box_muller.dat");

    for(int i = 0; i < (N_STEPS / 2); i++){
        x1 = prng(generator);
        x2 = prng(generator);
        x1 = sqrt(-log(1.0 - x1));
        x2 = 2 * pi * x2;
        y1 = AVERAGE + (x1 * cos(x2) * sqrt(2) * SIGMA);
        y2 = AVERAGE + (x1 * sin(x2) * sqrt(2) * SIGMA);

        file_b << y1 << endl << y2 << endl;
    }

    file_b.close();
}

double load_box_muller(int dim_subsample, vector<double>& y){
    /* Load Box-Muller data and compute the average */

    int steps = 0;
    double value, y_ave = 0.;

    ifstream data_bm("file_box_muller.dat");

    if (data_bm.is_open()){
        while ((data_bm >> value) && (steps <= dim_subsample)){
            steps++;
            y_ave += value;
            y.push_back(value);
        }
        y_ave = y_ave / steps;
        data_bm.close();
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }

    return y_ave;
}

inline double prf(double q, double q_try) {
    /* Probability ratio function for MC data */

    q = pow(q - AVERAGE, 2) - pow(q_try - AVERAGE, 2);
    return exp(q / (2*pow(SIGMA, 2)));
}

void generate_MC_data() {
    /* Generate MC data */

    double x, y, q_try, q = START, sample_average = 0.;

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

double load_and_average(int dim_therm, int dim_subsample, vector<double>& x){
    /* Load MC data and compute the average */

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

double blocking(int blocking_iteration, vector<double>& x){
    /* Compute sigma with the blocking algorithm */

    int steps = x.size(), k_steps = steps;
    double x_ave = 0., var = 0.;
    vector<double> x_k;

    x_k.reserve(steps);

    // compute average
    for(int i = 0; i < steps; i++) x_ave += x[i];
    x_ave = x_ave / steps;

    // compute variance
    for(int i = 0; i < steps; i++) var += pow(x[i] - x_ave, 2);
    var = var / (pow(steps, 2) - steps);

    // initialize working copy of x
    for(int i = 0; i < steps; i++) x_k.push_back(x[i]);

    // start blocking algorithm
    for(int k = 1; k <= blocking_iteration; k++){
        var = 0.;
        x_ave = 0.;
        k_steps = k_steps / 2;

        for(int i = 0; i < k_steps; i++){
            x_k[i] = (x_k[2*i] + x_k[2*i + 1]) / 2;
            x_ave += x_k[i];
        }
        x_ave = x_ave / k_steps;

        for(int i = 0; i < k_steps; i++) var += pow(x_k[i] - x_ave, 2);
        var = var / (pow(k_steps, 2) - k_steps);


        cout << "logging: k steps " << k << " -> ";
        cout << x_ave << " ± " << sqrt(var) << endl;
    }

  return sqrt(var);
}

double autocorrelation(vector<double>& x, float x_ave){
    /* Autocorrelation function, variance and error */

    int k = 0;
    int steps = x.size();
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

inline double estimator(vector<double>& x){
    /* Estimator to be evaluated on the sample*/

    int steps = x.size();
    double pow2 = 0., pow4 = 0.;

    for (int i = 0; i < steps; i++){
        pow2 += pow(x[i],2);
        pow4 += pow(x[i],4);
    }

    pow4 = (pow4 * steps) / (3 * pow2 * pow2);
    return pow4;
}

double bootstrap(vector<double>& x, int num_fake_samples){
    /* Bootstrap without autocorrelation */

    int steps = x.size();
    double value, fake_ave, sigma;
    vector<double> fake_sample, fake_values;

    fake_sample.reserve(steps);

    // Create the file for sigma(k)
    ofstream file_boot;
    file_boot.open("file_boot_w.dat");

    // Iterate over different k = number of fake samples
    for (int k = 10; k <= num_fake_samples; k += 10){

        fake_ave = 0.;
        fake_values.reserve(k);
        // Iterate over each fake sample
        for (int j = 0; j < k; j++){
            // Generate the j-fake_sample
            for (int i = 0; i < steps; i++){
                fake_sample.push_back(x[randomint(generator)]);
            }

            // Compute the estimator and store it
            value = estimator(fake_sample);
            fake_sample.clear();
            fake_ave += value;
            fake_values.push_back(value);
        }
        fake_ave = fake_ave / k;

        // Compute the variance of k fake samples
        sigma = 0.;
        for (auto val : fake_values) sigma += pow(val, 2);
        sigma = sigma / k;
        sigma = sigma - pow(fake_ave, 2);
        sigma = sqrt(sigma);

        fake_values.clear();

        cout << "logging: num fake samp " << k << " -> ";
        cout << fake_ave << " ± " << sigma << endl;

        // Save the sigma(k) to the file
        file_boot << k << " " << sigma << endl;
    }

    // Close and return sigma(num_fake_samples)
    file_boot.close();
    return sigma;
}

double bootstrap_corr(vector<double>& x, int num_fake_samples, int max_lenght){
    /* Bootstrap with autocorrelation */

    int steps = x.size(), draws, index, l;
    double value, fake_ave, sigma;
    vector<double> fake_sample, fake_values;

    fake_sample.reserve(steps);
    fake_values.reserve(num_fake_samples);

    // Create the file for sigma(k)
    ofstream file_boot;
    file_boot.open("file_boot_c.dat");

    // Iterate over different lenght of the correlated blocks
    for (int lenght = 1; lenght <= max_lenght; lenght *= 2){

        fake_ave = 0.;
        draws = 1 + ((steps - 1) / lenght);
        // Iterate over each fake sample
        for (int i = 0; i < num_fake_samples; i++){

            // Generate the i-fake_sample
            for (int j = 0; j < draws; j++){
                // Extract the j-random_int
                index = randomint(generator);
                // Push the j-block
                l = 0;
                while (l < lenght && fake_sample.size() < steps){
                    if (index + l == steps) index = index - steps;
                    fake_sample.push_back(x[index + l]);
                    l++;
                }
            }

            // Compute the estimator and store it
            value = estimator(fake_sample);
            fake_sample.clear();
            fake_ave += value;
            fake_values.push_back(value);
        }
        fake_ave = fake_ave / num_fake_samples;

        // Compute the variance of k fake samples
        sigma = 0.;
        for (auto val : fake_values) sigma += pow(val, 2);
        sigma = sigma / num_fake_samples;
        sigma = sigma - pow(fake_ave, 2);
        sigma = sqrt(sigma);

        fake_values.clear();

        cout << "logging: correl lenght " << lenght << " -> ";
        cout << fake_ave << " ± " << sigma << endl;

        // Save the sigma(lenght) to the file
        file_boot << lenght << " " << sigma << endl;
    }

    // Close and return sigma(circa max_lenght)
    file_boot.close();
    return sigma;
}

//---- Main --------------------------------------------------------------------

int main() {
    double x_ave, sigma;
    vector<double> data;

    // generate_BM_data();
    // generate_MC_data();

    x_ave = load_box_muller(DIM_SUBSAMP, data);
    // x_ave = load_and_average(DIM_THERM, DIM_SUBSAMP, data);
    cout << "The average of the sample is " << x_ave << endl;

    sigma = blocking(15, data);
    cout << "The result is " << x_ave << " +- " << sigma << endl;

    // sigma = autocorrelation(data, x_ave);
    // cout << "The result is " << x_ave << " +- " << sigma << endl;

    sigma = bootstrap(data, NUM_FAKESAMPLES);
    cout << "Bootstrap gives: " << sigma << endl;

    // sigma = bootstrap_corr(data, NUM_FAKESAMPLES, CORREL_LENGHT);
    // cout << "Bootstrap corr gives: " << sigma << endl;

    return 0;
}
