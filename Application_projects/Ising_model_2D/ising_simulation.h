/*******************************************************************************
*
* Subroutine for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#ifndef ISING_SIMULATION_H
#define ISING_SIMULATION_H

#include <iostream>
#include <fstream>
#include <string>

/* Import the Class lattice */
#include "ising_lattice.h"

/*
* CONFIGURATION PARAMETERS OF THE SIMULATION
* G_FLAG = lattice's geometry flag:
*          0 and others not implemented yet.
*          1 for 1D periodic chain.
*          2 for 2D square with PBC.
* I_FLAG = lattice's initial configuration flag:
*          0 for cold initialization.
*          1 for hot (random) initialization.
*          2 for loading the previous configuration.
* EXTFIELD = adimensional intensity of the external magnetic field.
* I_DECORREL = MC-updates of the lattice between different measurements.
* MEASURES = desired number of measures for each value of beta.
*/
#define G_FLAG 2
#define I_FLAG 1
#define EXTFIELD 0.
#define I_DECORREL 5
#define MEASURES 20000

using namespace std;

//----Contents------------------------------------------------------------------

/* MC-simulation for a given side and beta */
int run_simulation(int side, float beta){
    double ener, magn;
    string file_name;
    lattice ising(side, G_FLAG, I_FLAG);

    cout << "Start running with beta = " << beta << endl << endl;

    // Adapting the lattice to new beta
    for(int i = 0; i < 5 * I_DECORREL; i++) ising.update(beta, EXTFIELD);
    // Verify initial state
    magn = ising.magnetization();
    ener = ising.energy(EXTFIELD);

    // Create the file
    ofstream file;
    file_name = "Side_" + to_string(side) + "/side_" + to_string(side) +
                "_beta_" + to_string(beta) + ".dat";
    cout << "Creating file: " << file_name << endl << "-> starting energy = "
         << ener << endl << "-> starting magnet = " << magn << endl << endl;
    //ising.show_configuration();

    file.open(file_name);
    // Update ising and take measures
    for(int n = 0; n < MEASURES; n++){
        for(int i = 0; i < I_DECORREL; i++) ising.update(beta, EXTFIELD);

        ener = ising.energy(EXTFIELD);
        magn = ising.magnetization();

        file << ener << " " << magn << endl;
    }
    file.close();

    // We can stop the simulation when a file is completed
    //ising.save_configuration();
    cout << "Creation of " << file_name << " completed." << endl << endl;
    return 0;
}

#endif
