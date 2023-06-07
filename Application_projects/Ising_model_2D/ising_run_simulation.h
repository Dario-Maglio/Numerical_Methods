/*******************************************************************************
*
* Subroutine for the Ising simulation
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#ifndef ISING_RUN_SIMULATION_H
#define ISING_RUN_SIMULATION_H

#include <iostream>
#include <fstream>
#include <string>

// Import the Class lattice
#include "ising_lattice_class.h"

using namespace std;

/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* G_FLAG = lattice's geometry flag:
*          0 and others not implemented yet.
*          1 for 1D periodic chain.
*          2 for 2D square with PBC.
*
* I_FLAG = lattice's initial configuration flag:
*          0 for cold initialization.
*          1 for hot (random) initialization.
*          2 for loading the previous configuration and append data.
*
* EXTFIELD = adimensional intensity of the external magnetic field.
*
* I_DECORREL = MC-updates of the lattice between different measurements.
*
* MEASURES = desired number of measures for each value of beta.
*
*******************************************************************************/

#define G_FLAG 2
#define I_FLAG 2
#define EXTFIELD 0.
#define I_DECORREL 10
#define MEASURES 2000

//--- Contents -----------------------------------------------------------------

void run_simulation(int side, float beta){
    /* MC-simulation for a given side and beta */

    ofstream file;
    string directory, name_file, name_file_data, name_file_state, message;
    lattice ising(side, G_FLAG, I_FLAG);

    // Define path data directory
    directory = "Data_simulations/Side_" + to_string(side) + "/";
    // Define name file last configuration of the lattice
    name_file_state = "side_" + to_string(side) + "_beta_" + to_string(beta);
    // Define name file simulation for a given side and beta
    name_file_data =  name_file_state + ".dat";

    // Prepare the lattice for the simulation
    if(I_FLAG == 2){
        // Load last configuration of the lattice
        ising.load_configuration(directory + name_file_state);
    } else {
        // Thermalization phase
        for(int i = 0; i < (100*I_DECORREL); i++) ising.update(beta, EXTFIELD);
    }

    // Print initial energy and magnetization
    message = "File: " + name_file_data +
              "\n -> starting energy = " + to_string(ising.energy(EXTFIELD)) +
              "\n -> starting magnet = " + to_string(ising.magnetization());
    cout << message << endl << endl;

    // Open the data file
    if(I_FLAG == 2){
        file.open(directory + name_file_data, ios_base::app);
    } else {
        file.open(directory + name_file_data);
    }

    // Update ising and take measures
    for(int n = 0; n < MEASURES; n++){
        for(int i = 0; i < I_DECORREL; i++) ising.update(beta, EXTFIELD);
        file << ising.energy(EXTFIELD) << " " << ising.magnetization() << endl;
    }
    file.close();

    // We can stop the simulation when a file is completed
    ising.save_configuration(directory + name_file_state);
    cout << "Creation of " << name_file_data << " completed." << endl << endl;
}

#endif
