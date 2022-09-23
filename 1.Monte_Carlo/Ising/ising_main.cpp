/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>

/*
* import the the PRNG and the class lattice
*/
#include "ising_lattice.h"

/*
* SIDE = size of the lattice's side.
* G_FLAG = geometry flag; 1 for 1D periodic chain, 2 for 2D square with PBC.
* I_FLAG = initial configuration flag; 0 for cold initialization, 1 for hot
*   (random) initialization, 2 for loading the previous configuration from file.
*/
#define SIDE 25
#define G_FLAG 2
#define I_FLAG 1

using namespace std;

//----Contents------------------------------------------------------------------

int main(){
    // Note that we don't need to define random_number because we've imported it
    random_number = prng(generator);
    cout << random_number << endl << endl;

    lattice ising(SIDE, G_FLAG, I_FLAG);

    ising.show_configuration();

    ising.save_configuration();

    ising.show_nearest_neighbors();
}
