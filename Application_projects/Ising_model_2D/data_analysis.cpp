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

/*
* ALGORITHMS PARAMETERS
* BLOCKS_LENGHT = number of blocks reductions in the blocking algorithm
* CORR_LENGHT = lenght of the correlated blocks in the bootstrap algorithm
* NUM_FAKE_SAMP = number of fake samples of the bootstrap algorithm
*/
#define BLOCKS_LENGHT 7
#define MIN_CORR_LENGHT 25
#define MAX_CORR_LENGHT 200
#define NUM_FAKE_SAMP 300

// Define the namespace
using namespace std;

// Define the PRNG
#define SEED 42
mt19937 generator(SEED);
//random_device device;
//mt19937 generator(device());

//----Contents------------------------------------------------------------------

// Compute sigma with the blocking algorithm
double blocking(vector<double>& x, int blocking_iteration){
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
double bootstrap_corr(vector<double>& x, int num_fake_samples, int max_lenght){
    int steps = x.size(), lenght, draws, rand_index, k;
    double estimator, est_ave, est_var, samp_ave;
    vector<double> fake_sample, measures;

    fake_sample.reserve(steps);
    measures.reserve(num_fake_samples);
    uniform_int_distribution<long int> randomint(0, steps);

    // Iterate over different lenght of the correlated
    lenght = MIN_CORR_LENGHT;
    while(lenght <= max_lenght){

        est_ave = 0.;
        draws = 1 + ((steps - 1) / lenght);

        // Iterate over each fake sample
        for (int i = 0; i < num_fake_samples; i++){

            samp_ave = 0.;
            // generate the i-fake_sample
            for (int j = 0; j < draws; j++){
                // extract the j-random_int
                rand_index = randomint(generator);

                // push the j-block
                k = 0;
                while (k < lenght && fake_sample.size() < steps){
                    if (rand_index + k == steps) rand_index -= steps;
                    fake_sample.push_back(x[rand_index + k]);
                    samp_ave += x[rand_index + k];
                    k++;
                }
            }
            samp_ave = samp_ave / steps;

            // compute the variance, which here is our estimator
            estimator = 0.;
            for(auto val: fake_sample) estimator += pow(val - samp_ave, 2);
            estimator = estimator / (steps - 1);     // remember per volume

            fake_sample.clear();
            est_ave += estimator;
            measures.push_back(estimator);
        }
        est_ave = est_ave / num_fake_samples;

        // compute the variance over num_fake_samples fake samples
        est_var = 0.;
        for (auto val : measures) est_var += pow(val - est_ave, 2);
        est_var = est_var / (num_fake_samples - 1);  // per N divided by N - 1

        measures.clear();

        cout << "logging: correl lenght " << lenght << " -> ";
        cout << est_ave << " ± " << sqrt(est_var) << endl;
        lenght = (lenght * 4) / 3;
    }

    return sqrt(est_var);
}

/* Operations over each file */
void file_operations(int side, float beta, vector<double>& measures){
    int lenght = 0;
    double ene_ave = 0., mag_ave = 0., ene_var = 0., mag_var = 0., cumul = 0.;
    string file_name;
    vector<double> energies, magnetis;

    file_name = "Side_" + to_string(side) + "/side_" + to_string(side) +
                "_beta_" + to_string(beta) + ".dat";
    ifstream file(file_name);

    if (file.is_open()) {
        // Load data and compute averages
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

        // compute and store the averages
        ene_ave = ene_ave / lenght;
        mag_ave = mag_ave / lenght;
        measures[0] = ene_ave;
        measures[2] = mag_ave;

        // compute the errors with blocking algorithm
        cout << endl << "side: " << side << " | beta: " << beta << endl;
        cout << side << " " << beta << " --- energy error: " << endl;
        measures[1] = blocking(energies, BLOCKS_LENGHT);
        cout << side << " " << beta << " --- magnet error: " << endl;
        measures[3] = blocking(magnetis, BLOCKS_LENGHT);

        // compute the estimators
        for(auto val: energies) ene_var += pow(val - ene_ave, 2);
        measures[4] = ene_var / (lenght - 1);    // remember per volume

        for(auto val: magnetis) mag_var += pow(val - mag_ave, 2);
        measures[6] = mag_var / (lenght - 1);    // remember per volume

        // compute the Binder cumulant
        mag_var = 0.;
        for(auto val: magnetis) mag_var += pow(val, 2);
        for(auto val: magnetis) cumul += pow(val, 4);
        cumul = (cumul * lenght) / pow(mag_var, 2);
        measures[8] = cumul;

        // compute the error with bootstrap algorithm
        cout << side << " " << beta << " --- heat error: " << endl;
        measures[5] = bootstrap_corr(energies, NUM_FAKE_SAMP, MAX_CORR_LENGHT);
        cout << side << " " << beta << " --- chi error: " << endl;
        measures[7] = bootstrap_corr(magnetis, NUM_FAKE_SAMP, MAX_CORR_LENGHT);
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
}

void complete_analysis(){
    ofstream file;
    string file_name;
    vector<double> measures;

    measures.reserve(9);
    for(int i = 0; i < 9; i++) measures.push_back(0.);

    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){

        file_name = "Side_" + to_string(side) + "/averages_and_variance.dat";
        cout << endl << endl << "Creating file: " << file_name << endl;
        file.open(file_name);

        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta, measures);
            file << beta << " ";
            for(int i = 0; i < 9; i++) file << measures[i] << " ";
            file << endl;
        }
        file.close();
    }
}

void partial_analysis(){
    int side_in;
    ofstream file;
    string file_name;
    vector<double> measures;

    measures.reserve(9);
    for(int i = 0; i < 9; i++) measures.push_back(0.);
    
    cout << "Insert side lenght: ";
    cin >>  side_in; 

    file_name = "Side_" + to_string(side_in) + "/averages_and_variance.dat";
    cout << endl << "Creating file: " << file_name << endl;
    file.open(file_name);

    for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
        file_operations(side_in, beta, measures);
        file << beta << " ";
        for(int i = 0; i < 9; i++) file << measures[i] << " ";
        file << endl;
    }
    file.close();
}

/* Main for the data analysis */
int main(){
    
    // complete_analysis();

    partial_analysis();

    cout << "The work is done." << endl << endl;
}
