/*******************************************************************************
*
* Main program for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <fstream>
#include <random>        // import mt19937 periodo di 2^19937 − 1
#include <cmath>

#include "ising_lattice.h"

/*
* SIDE = size of the lattice's side.
* G_FLAG = geometry flag; 0 for 1D periodic chain, 1 for 2D square with PBC.
* I_FLAG = initial configuration flag; 0 for cold initialization, 1 for hot
*   (random) initialization, 2 for loading the previous configuration from file.
*/
#define SIDE 5
#define G_FLAG 1
#define I_FLAG 1

using namespace std;

//----Contents------------------------------------------------------------------

int main(){
    lattice ising(SIDE, G_FLAG, I_FLAG);

    ising.show_configuration();

    ising.save_configuration();

    // ising.show_nearest_neighbors();
}
