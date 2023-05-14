/*******************************************************************************
*
* Data analysis program for the outcomes of the Ising simulations
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <random>

using namespace std;

// Define the PRNG
#define SEED 42
mt19937 generator(SEED);
//random_device device;
//mt19937 generator(device());

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* BETA_SEP = separation between the betas of different simulations.
*
* SIDE_SEP = separation between the sides of different simulations.
*
*******************************************************************************/

#define DATA 9
#define SIDE_MIN 20
#define SIDE_MAX 60
#define SIDE_SEP 10
#define BETA_INI 0.3600
#define BETA_FIN 0.5100
#define BETA_SEP 0.0025

/*******************************************************************************
* PARAMETERS OF THE ALGORITHMS
*
* BLOCKING = number of blocks reductions in the blocking algorithm.
*
* CORR_LENGHT = lenght of the correlated blocks in the bootstrap algorithm.
*
* NUM_FAKE_SAMP = number of fake samples of the bootstrap algorithm.
*
*******************************************************************************/

#define BLOCKING 7
#define MIN_CORR_LENGHT 2
#define MAX_CORR_LENGHT 300
#define NUM_FAKE_SAMP 300

//--- Contents -----------------------------------------------------------------

double blocking(vector<double>& x, int blocking_iteration, ofstream &file){
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


        file << "k steps " << k << " -> ";
        file << x_ave << " ± " << sqrt(var) << endl;
    }

  return sqrt(var);
}

double bootstrap_corr(vector<double>& x, int num_fake_samples,
                      int max_lenght, ofstream &file){
    /* Bootstrap with autocorrelation */

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
            // Generate the i-fake_sample
            for (int j = 0; j < draws; j++){
                // Extract the j-random_int
                rand_index = randomint(generator);

                // Push the j-block
                k = 0;
                while (k < lenght && fake_sample.size() < steps){
                    if (rand_index + k == steps) rand_index -= steps;
                    fake_sample.push_back(x[rand_index + k]);
                    samp_ave += x[rand_index + k];
                    k++;
                }
            }
            samp_ave = samp_ave / steps;

            // Compute the variance, which here is our estimator
            estimator = 0.;
            for(auto val: fake_sample) estimator += pow(val - samp_ave, 2);
            estimator = estimator / (steps - 1);     // remember per volume

            fake_sample.clear();
            est_ave += estimator;
            measures.push_back(estimator);
        }
        est_ave = est_ave / num_fake_samples;

        // Compute the variance over num_fake_samples fake samples
        est_var = 0.;
        for (auto val : measures) est_var += pow(val - est_ave, 2);
        est_var = est_var / (num_fake_samples - 1);  // per N divided by N - 1

        measures.clear();

        file << "correl lenght " << lenght << " -> ";
        file << est_ave << " ± " << sqrt(est_var) << endl;
        lenght = lenght * 2;
    }

    return sqrt(est_var);
}

void file_operations(int side, float beta, vector<double>& measures,
                     ofstream &file_analysis){
    /* Operations over each file */

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

        // Compute and store the averages
        ene_ave = ene_ave / lenght;
        mag_ave = mag_ave / lenght;
        measures[0] = ene_ave;
        measures[2] = mag_ave;

        // Compute the errors with blocking algorithm
        cout << endl << "side: " << side << " | beta: " << beta << endl;
        file_analysis << endl << "L: " << side << " | beta: " << beta << endl;
        file_analysis << side << " " << beta << " --- energy error: " << endl;
        measures[1] = blocking(energies, BLOCKING, file_analysis);
        file_analysis << side << " " << beta << " --- magnet error: " << endl;
        measures[3] = blocking(magnetis, BLOCKING, file_analysis);

        // Compute the estimators
        for(auto val: energies) ene_var += pow(val - ene_ave, 2);
        measures[4] = ene_var / (lenght - 1);    // remember per volume

        for(auto val: magnetis) mag_var += pow(val - mag_ave, 2);
        measures[6] = mag_var / (lenght - 1);    // remember per volume

        // Compute the Binder cumulant
        mag_var = 0.;
        for(auto val: magnetis) mag_var += pow(val, 2);
        for(auto val: magnetis) cumul += pow(val, 4);
        cumul = (cumul * lenght) / pow(mag_var, 2);
        measures[8] = cumul;

        // Compute the error with bootstrap algorithm
        file_analysis << side << " " << beta << " --- heat error: " << endl;
        measures[5] = bootstrap_corr(energies, NUM_FAKE_SAMP, MAX_CORR_LENGHT,
                                     file_analysis);
        file_analysis << side << " " << beta << " --- chi error: " << endl;
        measures[7] = bootstrap_corr(magnetis, NUM_FAKE_SAMP, MAX_CORR_LENGHT,
                                     file_analysis);
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
}

void complete_analysis(){
    string file_name;
    ofstream file, file_analysis;
    vector<double> measures;

    measures.reserve(DATA);
    for(int i = 0; i < DATA; i++) measures.push_back(0.);

    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        file_name = "Data_analysis/side_" + to_string(side) + "_analysis.txt";
        cout << endl << "Creating file: " << file_name << endl;
        file_analysis.open(file_name);

        file_name = "Data_analysis/side_" + to_string(side) + "_data.dat";
        cout << endl << "Creating file: " << file_name << endl;
        file.open(file_name);

        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta, measures, file_analysis);
            file << beta << " ";
            for(int i = 0; i < DATA; i++) file << measures[i] << " ";
            file << endl;
        }
        file.close();
        file_analysis.close();
    }
}

void partial_analysis(){
    int side_in;
    string file_name;
    ofstream file, file_analysis;
    vector<double> measures;

    measures.reserve(DATA);
    for(int i = 0; i < DATA; i++) measures.push_back(0.);

    cout << "Insert side lenght: ";
    cin >>  side_in;

    file_name = "Data_analysis/side_" + to_string(side_in) + "_analysis.txt";
    cout << endl << "Creating file: " << file_name << endl;
    file_analysis.open(file_name);

    file_name = "Data_analysis/side_" + to_string(side_in) + "_data.dat";
    cout << endl << "Creating file: " << file_name << endl;
    file.open(file_name);

    for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
        file_operations(side_in, beta, measures, file_analysis);
        file << beta << " ";
        for(int i = 0; i < DATA; i++) file << measures[i] << " ";
        file << endl;
    }
    file.close();
    file_analysis.close();
}

//--- Main ---------------------------------------------------------------------

int main(){
    /* Main program for the data analysis. */

    // complete_analysis();

    partial_analysis();

    cout << "The work is done." << endl << endl;
}
