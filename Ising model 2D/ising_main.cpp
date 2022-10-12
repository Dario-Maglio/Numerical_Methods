/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <fstream>
#include <string>

/* Import the Class lattice */
#include "ising_lattice.h"

/*
* CONFIGURATION PARAMETERS OF THE LATTICE
* SIDE = size of the lattice's side.
* G_FLAG = geometry flag; 1 for 1D periodic chain, 2 for 2D square with PBC.
* I_FLAG = initial configuration flag; 0 for cold initialization, 1 for hot
*   (random) initialization, 2 for loading the previous configuration from file.
*/
#define SIDE 20
#define G_FLAG 2
#define I_FLAG 1

/*
* CONFIGURATION PARAMETERS OF THE PHYSICS OF THE PROBLEM
* BETA = reciprocal of the product of the temperature and k_B.
* EXTFIELD = adimensional intensity of the external magnetic field.
*/
#define INITIAL_BETA 0.345
#define FINAL_BETA 0.555
#define EXTFIELD 0.

/*
* CONFIGURATION PARAMETERS OF THE MEASUREMENTS
* I_DECORREL = updates of the system between different measurements.
* MEASURES = desired number of measures.
*/
#define I_DECORREL 10
#define MEASURES 10000

using namespace std;

//----Contents------------------------------------------------------------------

int main(){
    float beta;
    double ener, magn;
    string file_name;
    lattice ising(SIDE, G_FLAG, I_FLAG);

    for(float beta = INITIAL_BETA; beta <= FINAL_BETA; beta += 0.01){
        cout << "Running with beta = " << beta << endl;

        // Adapting the lattice to new beta
        for(int i = 0; i < I_DECORREL; i++) ising.update(beta, EXTFIELD);

        ener = ising.energy(EXTFIELD);
        cout << "-> starting energy = " << ener << endl;
        magn = ising.magnetization();
        cout << "-> starting magnet = " << magn << endl;
        ising.show_configuration();

        // Create the file
        ofstream file;
        file_name = "Side_" + to_string(SIDE) + "/side_" + to_string(SIDE) +
                    "_beta_" + to_string(beta) + ".dat";
        cout << endl << "Creating file: " << file_name << endl;
        file.open(file_name);

        // Update ising and take measures
        for(int n = 0; n < MEASURES; n++){
            for(int i = 0; i < I_DECORREL; i++) ising.update(beta, EXTFIELD);

            ener = ising.energy(EXTFIELD);
            magn = ising.magnetization();

            file << ener << " " << magn << endl;
        }

        file.close();
        cout << "The work is done." << endl << endl;
    }

    ising.save_configuration();
}
