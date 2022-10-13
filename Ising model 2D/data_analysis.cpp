/*******************************************************************************
*
* Data analysis program for the outcomes of the Ising simulations
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <fstream>
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
void file_operations(int side, float beta){
    string file_name;
    file_name = "Side_" + to_string(side) + "/side_" + to_string(side) +
                "_beta_" + to_string(beta) + ".dat";
    cout << endl << "Loading file: " << file_name << endl;

    ifstream file(file_name);
    if (file.is_open()) {
        double ene, mag;
        while (file >> ene >> mag){
            cout << "--> energy " << ene << " | magnet " << mag << endl;;
        }
        file.close();
    } else {
        cerr << "Error: unable to open the file." << endl;
        exit(1);
    }
}

/* Main for the data analysis */
int main(){

    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            file_operations(side, beta);
        }
    }

    cout << "The work is done." << endl << endl;

}
