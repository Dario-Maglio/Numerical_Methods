/**
* Main program for the Ising simulation bla bla bla
*/

//----Preprocessor directives--------------------------------------
#include <iostream>
#include <fstream>       // file stream
#include <random>        // import mt19937 periodo di 2^19937 − 1
#include <cmath>

#include "ising_lattice.h"

#define SEED 42

#define L 5
#define G_FLAG 0
#define T_FLAG 0

using namespace std;

//----Contents-----------------------------------------------------

int main(){
    lattice ising(L, G_FLAG, T_FLAG);

    ising.show_configuration();

    ising.show_nearest_neighbors();
}
