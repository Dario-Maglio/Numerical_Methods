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
mt19937_64 generator(SEED);
//random_device device;
//mt19937_64 generator(device());

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* BETA_SEP = separation between the betas of different simulations.
*
* SIDE_SEP = separation between the sides of different simulations.
*
*******************************************************************************/

#define DATA 10
#define SIDE_MIN 20
#define SIDE_MAX 60
#define SIDE_SEP 10

#define BETA_INI 0.3600
#define BETA_FIN 0.5100
#define BETA_SEP 0.0025
#define BETA_C_INI 0.4153
#define BETA_C_FIN 0.4500
#define BETA_C_SEP 0.0010

/*******************************************************************************
* PARAMETERS OF THE ANALYSIS
*
* BLOCKS = number of block-reductions in the blocking algorithm.
*
* MIN_CORR_LENGHT = min dim of the correlated blocks in the bootstrap algorithm.
*
* MAX_CORR_LENGHT = max dim of the correlated blocks in the bootstrap algorithm.
*
* NUM_FAKE_SAMP = number of fake samples in the bootstrap algorithm.
*
* DIM_FAKE_SAMP = dimension of the fake samples in the bootstrap algorithm.
*
*******************************************************************************/

#define BLOCKS 6
#define MIN_CORR_LENGHT 2
#define MAX_CORR_LENGHT 512
#define NUM_FAKE_SAMP 300
#define DIM_FAKE_SAMP 150000

//--- Contents -----------------------------------------------------------------

double estimator(vector<double>& x){
    /* Susceptibility and specific heat / volume */

    int dim = x.size();
    double average = 0., estimator = 0;

    for(auto val: x) average += val;
    average = average / dim;

    for(auto val: x) estimator += pow(val - average, 2);
    estimator = estimator / (dim - 1);
    return estimator;
}

double binder(vector<double>& x){
    /* Binder cumulant function */

    double x_pow2 = 0., x_pow4 = 0.;

    for(auto val: x) x_pow4 += pow(val, 4);
    for(auto val: x) x_pow2 += pow(val, 2);
    return ((x_pow4 * x.size())/ pow(x_pow2, 2));
}

double blocking(vector<double>& x, ofstream &file){
    /* Compute sigma with the blocking algorithm */
    // Given a sample x, compute its average and variance while varying the
    // block's dimension from 1 to 2^BLOCKS and writing them
    // into a file. The function blocking returns the last std deviation.

    int steps = x.size(), k_steps;
    double x_ave, var;
    vector<double> x_k;

    // initialize working copy of x
    x_k.reserve(steps);
    for(int i = 0; i < steps; i++) x_k.push_back(x[i]);
    file << "Initial sample dimension: " << steps << endl;

    // start blocking algorithm
    for(int k = 1; k <= BLOCKS; k++){
        var = 0.;
        x_ave = 0.;
        k_steps = steps / pow(2, k);

        for(int i = 0; i < k_steps; i++){
            x_k[i] = (x_k[2*i] + x_k[2*i + 1]) / 2;
            x_ave += x_k[i];
        }
        x_ave = x_ave / k_steps;

        // compute variance
        for(int i = 0; i < k_steps; i++) var += pow(x_k[i] - x_ave, 2);
        var = var / (k_steps * (k_steps - 1));

        file << "k steps " << k << " -> ";
        file << x_ave << " ± " << sqrt(var) << endl;
    }

    file << "Final sample dimension: " << k_steps << endl;
    return sqrt(var);
}

double bootstrap_corr(vector<double>& x, double (*estimator)(vector<double>& ), ofstream &file){
    /* Compute sigma with the bootstrap algorithm */
    // Given a sample x, generate NUM_FAKE_SAMP samples via resampling
    // and compute average and variance of an estimator while varying the
    // correlation lenght from MIN to MAX_CORR_LENGHT and writing them
    // into a file. The function bootstrap returns the last std deviation.

    int rand_index;
    double est_val, est_ave, est_var;
    vector<double> fake_sample, measures;

    measures.reserve(NUM_FAKE_SAMP);
    fake_sample.reserve(DIM_FAKE_SAMP);
    uniform_int_distribution<long int> randomint(0, x.size());

    // iterate over different correlation lenghts
    for(int lenght = MIN_CORR_LENGHT; lenght <= MAX_CORR_LENGHT; lenght *= 2){

        // initialite the average of the
        // estimator over the fake samples
        est_ave = 0.;

        // evaluate the estimator over each fake sample
        // and save the result into the vector measures
        for(int i = 0; i < NUM_FAKE_SAMP; i++){

            // generate the i-fake_sample
            for (int j = 0; fake_sample.size() < DIM_FAKE_SAMP; j++){
                // extract the j-random_int and push the j-block
                rand_index = randomint(generator);
                for (int k = 0; k < lenght; k++){
                    if (rand_index + k == x.size()) break;
                    if (fake_sample.size() == DIM_FAKE_SAMP) break;
                    fake_sample.push_back(x[rand_index + k]);
                }
            }

            // evaluate the estimator and save the value
            est_val = estimator(fake_sample);
            est_ave += est_val;
            measures.push_back(est_val);
            fake_sample.clear();
        }
        est_ave = est_ave / NUM_FAKE_SAMP;

        // compute the estimator variance over fake samples
        est_var = 0.;
        for (auto val : measures) est_var += pow(val - est_ave, 2);
        est_var = est_var / (NUM_FAKE_SAMP - 1);

        measures.clear();

        file << "correl lenght " << lenght << " -> ";
        file << est_ave << " ± " << sqrt(est_var) << endl;
    }

    return sqrt(est_var);
}

//--- File operations ----------------------------------------------------------

void file_operations(int side, float beta, vector<double>& data, ofstream &file_analysis){
    /* Operations over each beta file */
    // For a given side and beta, compute average energy, average
    // magnetization, susceptibility, specific heat and errors.

    int measures = 0;
    double ene, mag, ene_ave = 0., mag_ave = 0.;
    string file_name;
    vector<double> energies, magnetis;

    file_name = "Data_simulations/Side_" + to_string(side) + "/side_";
    file_name += to_string(side) + "_beta_" + to_string(beta) + ".dat";
    ifstream file(file_name);

    if (file.is_open()) {
        // load data and compute averages
        while (file >> ene >> mag){
            measures++;
            mag = abs(mag);
            ene_ave += ene;
            mag_ave += mag;
            energies.push_back(ene);
            magnetis.push_back(mag);
        }
        file.close();

        cout << "side: " << side << " | beta: " << beta << endl << endl;
        file_analysis << "-----------------------------" << endl;
        file_analysis << "side: " << side << " | beta: " << beta << endl;
        file_analysis << "-----------------------------" << endl;

        // compute and store the averages
        ene_ave = ene_ave / measures;
        mag_ave = mag_ave / measures;
        data[0] = ene_ave;
        data[2] = mag_ave;

        // compute the errors with blocking algorithm
        file_analysis << " --- energy error: " << endl;
        data[1] = blocking(energies, file_analysis);
        file_analysis << " --- magnet error: " << endl;
        data[3] = blocking(magnetis, file_analysis);

        // Compute specific heat and susceptibility
        data[4] = estimator(energies) * pow(side, 2);
        data[6] = estimator(magnetis) * pow(side, 2);

        // compute the errors with bootstrap algorithm
        file_analysis << " --- heat error / V : " << endl;
        data[5] = bootstrap_corr(energies, &estimator, file_analysis);
        data[5] = data[5] * pow(side, 2);
        file_analysis << " --- chi error / V : " << endl;
        data[7] = bootstrap_corr(magnetis, &estimator, file_analysis);
        data[7] = data[7] * pow(side, 2);

        // compute the Binder cumulant
        file_analysis << " --- cumulant error: " << endl;
        data[8] = binder(magnetis);
        data[9] = bootstrap_corr(magnetis, &binder, file_analysis);

        file_analysis << "-----------------------------" << endl << endl;
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
}

void partial_analysis(){
    /* Analysis for a specific side */
    // Given an input side from bash, call file
    // operations for each beta selected above.

    int side_in;
    string file_name;
    ofstream file_data, file_analysis;
    vector<double> data;

    data.reserve(DATA);
    for(int i = 0; i < DATA; i++) data.push_back(0.);

    cout << "Insert side lenght: ";
    cin >>  side_in;
    cout << endl;

    file_name = "Data_analysis/side_" + to_string(side_in) + "_analysis.txt";
    cout << "Creating file: " << file_name << endl;
    file_analysis.open(file_name);

    file_name = "Data_analysis/side_" + to_string(side_in) + "_data.dat";
    cout << endl << "Creating file: " << file_name << endl << endl;
    file_data.open(file_name);

    // begin loop over betas
    for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
        file_operations(side_in, beta, data, file_analysis);
        file_data << beta << " ";
        for(int i = 0; i < DATA; i++) file_data << data[i] << " ";
        file_data << endl;
    }
    // begin loop over central betas
    file_analysis << "\n\nBegin loop over central betas\n";
    for(float beta = BETA_C_INI; beta <= BETA_C_FIN; beta += BETA_C_SEP){
        file_operations(side_in, beta, data, file_analysis);
        file_data << beta << " ";
        for(int i = 0; i < DATA; i++) file_data << data[i] << " ";
        file_data << endl;
    }

    file_data.close();
    file_analysis.close();
}

void complete_analysis(){
    /* Analysis for each side */
    // Call file operations for each
    // side and beta selected above.

    string file_name;
    ofstream file_data, file_analysis;
    vector<double> data;

    data.reserve(DATA);
    for(int i = 0; i < DATA; i++) data.push_back(0.);

    // begin loop over sides
    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        file_name = "Data_analysis/side_" + to_string(side) + "_analysis.txt";
        cout << endl << "Creating file: " << file_name << endl;
        file_analysis.open(file_name);

        file_name = "Data_analysis/side_" + to_string(side) + "_data.dat";
        cout << endl << "Creating file: " << file_name << endl << endl;
        file_data.open(file_name);

        // begin loop over betas
        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta, data, file_analysis);
            file_data << beta << " ";
            for(int i = 0; i < DATA; i++) file_data << data[i] << " ";
            file_data << endl;
        }
        // begin loop over central betas
        file_analysis << "\n\nBegin loop over central betas\n";
        for(float beta = BETA_C_INI; beta <= BETA_C_FIN; beta += BETA_C_SEP){
            file_operations(side, beta, data, file_analysis);
            file_data << beta << " ";
            for(int i = 0; i < DATA; i++) file_data << data[i] << " ";
            file_data << endl;
        }

        file_data.close();
        file_analysis.close();
    }
}

//--- Main ---------------------------------------------------------------------

int main(){
    /* Main program for the data analysis. */

    //partial_analysis();
    complete_analysis();

    cout << "The work is done." << endl << endl;
}
