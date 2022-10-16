/*******************************************************************************
*
* Data analysis program for the outcomes of the Ising simulations
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <random>

/*
* CONFIGURATION PARAMETERS
* SIDE_SEP = separation between the sides of different simulations.
* BETA_SEP = separation between the betas of different simulations.
*/
#define SIDE_MIN 20
#define SIDE_MAX 60
#define SIDE_SEP 10
#define BETA_INI 0.3600
#define BETA_FIN 0.5100
#define BETA_SEP 0.0025
#define SEED 42
#define BLOCKS 5
#define NUM_FAKE_SAMPLES 100
#define MAX_CORR_LENGHT 200

using namespace std;

// define the PRNG
mt19937_64 generator(SEED);

//----Contents------------------------------------------------------------------

// Compute sigma with the blocking algorithm
double blocking(int blocking_iteration, vector<double>& x){
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

/* Bootstrap with autocorrelation */
double bootstrap_corr(int num_fake_samples, int max_lenght, vector<double>& x){
    int steps = x.size(), draws, index, l;
    double estimator, fake_ave, sigma, ave;
    vector<double> fake_sample, fake_values;

    fake_sample.reserve(steps);
    fake_values.reserve(num_fake_samples);
    uniform_int_distribution<long int> randomint(0, steps);

    // Iterate over different lenght of the correlated blocks
    for (int lenght = 8; lenght <= max_lenght; lenght *= 2){

        fake_ave = 0.;
        draws = 1 + ((steps - 1) / lenght);
        // Iterate over each fake sample
        for (int i = 0; i < num_fake_samples; i++){
  
            estimator = 0.;
            ave = 0.;
            // Generate the i-fake_sample
            for (int j = 0; j < draws; j++){
                // Extract the j-random_int
                index = randomint(generator);

                // Push the j-block
                l = 0;
                while (l < lenght && fake_sample.size() < steps){
                    if (index + l == steps) index = index - steps;
                    fake_sample.push_back(x[index + l]);
                    ave += x[index + l]; 
                    l++;
                }
            }

            // compute variance
            ave = ave / steps;
            for(int i = 0; i < steps; i++) estimator += pow(fake_sample[i] - ave, 2);
            estimator = estimator / (pow(steps, 2) - steps);

            fake_sample.clear();
            fake_ave += estimator;
            fake_values.push_back(estimator);
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

    }


    return sigma;
}

/* Operations over each file */
void file_operations(int side, float beta, vector<double>& measures){
    double ene_ave = 0., mag_ave = 0., ene_var = 0., mag_var = 0., cumul = 0.;
    string file_name;
    vector<double> energies, magnetis;

    file_name = "Side_" + to_string(side) + "/side_" + to_string(side) +
                "_beta_" + to_string(beta) + ".dat";
    ifstream file(file_name);

    if (file.is_open()) {
        int lenght = 0;
        double ene, mag;

        while (file >> ene >> mag){
            lenght++;
            mag = abs(mag);
            ene_ave += ene;
            mag_ave += mag;
            energies.push_back(ene);
            magnetis.push_back(mag);
        }

        file.close();

        // compute the averages
        ene_ave = ene_ave / lenght;
        mag_ave = mag_ave / lenght;
        measures[0] = ene_ave;
        measures[2] = mag_ave;

        // compute the error with blocking algorithm
        cout << side << " " << beta << " ---energy error: " << endl;
        measures[1] = blocking(BLOCKS, energies);
        cout << side << " " << beta << " ---magnetization error: " << endl;
        measures[3] = blocking(BLOCKS, magnetis);

        // compute the variance and Binder cumulant
        for(auto val: energies) ene_var += pow(val, 2);
        ene_var = (ene_var / lenght) - pow(ene_ave, 2);
        measures[4] = ene_var;

        for(auto val: magnetis) mag_var += pow(val, 2);
        for(auto val: magnetis) cumul += pow(val, 4);
        cumul = (cumul * lenght) / pow(mag_var, 2);
        measures[8] = cumul;

        mag_var = mag_var / lenght - pow(mag_ave, 2);
        measures[6] = mag_var;

        // compute the error with bootstrap algorithm 
        cout << side << " " << beta << " --- C error: " << endl;
        measures[5] = bootstrap_corr(NUM_FAKE_SAMPLES, MAX_CORR_LENGHT, energies);
        cout << side << " " << beta << " --- Ki error: " << endl;
        measures[7] = bootstrap_corr(NUM_FAKE_SAMPLES, MAX_CORR_LENGHT, magnetis);

    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
}

/* Main for the data analysis */
int main(){
    string file_name;
    vector<double> measures;
    ofstream file;

    measures.reserve(9);
    for(int i = 0; i < 9; i++) measures.push_back(0.);

    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){

        file_name = "Side_" + to_string(side) + "/averages_and_variance.dat";
        cout << endl << "Creating file: " << file_name << endl;
        file.open(file_name);

        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta, measures);
            file << beta << " ";
            for(int i = 0; i < 9; i++) file << measures[i] << " ";
            file << endl;
        }

        file.close();
    }

    cout << "The work is done." << endl << endl;

}
