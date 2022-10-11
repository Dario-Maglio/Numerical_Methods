/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <fstream>

/*
* Import the class lattice
*/
#include "ising_lattice.h"

/*
* CONFIGURATION PARAMETERS OF THE LATTICE
* SIDE = size of the lattice's side.
* G_FLAG = geometry flag; 1 for 1D periodic chain, 2 for 2D square with PBC.
* I_FLAG = initial configuration flag; 0 for cold initialization, 1 for hot
*   (random) initialization, 2 for loading the previous configuration from file.
*/
#define SIDE 10
#define G_FLAG 2
#define I_FLAG 1

/*
* CONFIGURATION PARAMETERS OF THE PHYSICS OF THE PROBLEM
* BETA = reciprocal of the product of the temperature and k_B.
* EXTFIELD = adimensional intensity of the external magnetic field.
*/
#define BETA 0.02
#define EXTFIELD 0.

/*
* CONFIGURATION PARAMETERS OF THE MEASUREMENTS
* I_DECORREL = updates of the system between different measurements.
* MEASURES = desired number of measures.
*/
#define I_DECORREL 1000
#define MEASURES 1000

using namespace std;

//----Contents------------------------------------------------------------------

int main(){
    double ener, magn;
    lattice ising(SIDE, G_FLAG, I_FLAG);

    ofstream file;
    file.open("lattice_measures.dat");

    ener = ising.energy(EXTFIELD);
    cout << "The energy is " << ener << endl;
    magn = ising.magnetization();
    cout << "The magnetization is " << magn << endl;
    ising.show_configuration();

    for(int n = 0; n < MEASURES; n++){
        for(int i = 0; i < I_DECORREL; i++) ising.update(BETA, EXTFIELD);

        ener = ising.energy(EXTFIELD);
        magn = ising.magnetization();

        file << ener << " " << magn << endl;
    }

    file.close();

    ener = ising.energy(EXTFIELD);
    cout << "The energy is " << ener << endl;
    magn = ising.magnetization();
    cout << "The magnetization is " << magn << endl;
    ising.show_configuration();

    ising.save_configuration();
}
