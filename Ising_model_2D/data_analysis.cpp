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

using namespace std;

//----Contents------------------------------------------------------------------

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
        measures[1] = mag_ave;

        // compute the variance and Binder cumulant
        for(auto val: energies) ene_var += pow(val, 2);
        ene_var = (ene_var / lenght) - pow(ene_ave, 2);
        measures[2] = ene_var;

        for(auto val: magnetis) mag_var += pow(val, 2);
        for(auto val: magnetis) cumul += pow(val, 4);
        cumul = (cumul * lenght) / pow(mag_var, 2);
        measures[4] = cumul;

        mag_var = mag_var / lenght - pow(mag_ave, 2);
        measures[3] = mag_var;

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

    measures.reserve(5);
    for(int i = 0; i < 5; i++) measures.push_back(0.);

    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){

        file_name = "Side_" + to_string(side) + "/averages_and_variance.dat";
        cout << endl << "Creating file: " << file_name << endl;
        file.open(file_name);

        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta, measures);
            file << beta << " ";
            for(int i = 0; i < 5; i++) file << measures[i] << " ";
            file << endl;
        }

        file.close();
    }

    cout << "The work is done." << endl << endl;

}
