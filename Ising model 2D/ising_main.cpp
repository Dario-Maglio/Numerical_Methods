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
* SIDE_SEP = separation between the sides of different simulations.
* G_FLAG = lattice's geometry flag:
*          0 and others not implemented yet.
*          1 for 1D periodic chain.
*          2 for 2D square with PBC.
* I_FLAG = lattice's initial configuration flag:
*          0 for cold initialization.
*          1 for hot (random) initialization.
*          2 for loading the previous configuration.
*/
#define SIDE_MIN 20
#define SIDE_MAX 60
#define SIDE_SEP 10
#define G_FLAG 2
#define I_FLAG 1

/*
* CONFIGURATION PARAMETERS OF THE SIMULATION
* BETA_SEP = separation between the betas of different simulations.
* EXTFIELD = adimensional intensity of the external magnetic field.
* I_DECORREL = MC-updates of the lattice between different measurements.
* MEASURES = desired number of measures for each value of beta.
*/
#define BETA_INI 0.3600
#define BETA_FIN 0.5100
#define BETA_SEP 0.0025
#define EXTFIELD 0.
#define I_DECORREL 5
#define MEASURES 1000

using namespace std;

//----Contents------------------------------------------------------------------

/* MC-simulation for a given lattice and beta */
void run_simulation(lattice &ising, float beta){
    int side = ising.side_lenght;
    double ener, magn;
    string file_name;

    cout << "Running with beta = " << beta << endl;

    // Adapting the lattice to new beta
    for(int i = 0; i < I_DECORREL; i++) ising.update(beta, EXTFIELD);

    ener = ising.energy(EXTFIELD);
    cout << "-> starting energy = " << ener << endl;
    magn = ising.magnetization();
    cout << "-> starting magnet = " << magn << endl << endl;
    // ising.show_configuration();

    // Create the file
    ofstream file;
    file_name = "Side_" + to_string(side) + "/side_" + to_string(side) +
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

    // Close the file and save configuration
    // We can stop the simulation when a file is completed
    file.close();
    ising.save_configuration();
    cout << "Creation completed." << endl << endl;
}

/* Main program iterates over sides and betas */
int main(){
    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        lattice ising(side, G_FLAG, I_FLAG);

        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            run_simulation(ising, beta);
        }
    }
    cout << "The work is done." << endl << endl;
}
