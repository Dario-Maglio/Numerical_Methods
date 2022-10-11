/*******************************************************************************
*
* Test program for the Ising simulation
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>

/*
* import the class lattice
*/
#include "ising_lattice.h"

/*
* SIDE = size of the lattice's side.
* G_FLAG = geometry flag; 1 for 1D periodic chain, 2 for 2D square with PBC.
* I_FLAG = initial configuration flag; 0 for cold initialization, 1 for hot
*   (random) initialization, 2 for loading the previous configuration from file.
*/
#define SIDE 10
#define G_FLAG 2
#define I_FLAG 1

/*
* BETA = reciprocal of the temperature per k_B.
* EXTFIELD = intensity of the external magnetic field.
*/
#define BETA 0.0002
#define EXTFIELD 0.

#define SUBSAMP 10

using namespace std;

//----Contents------------------------------------------------------------------

int main(){
    double ener, magn;
    lattice ising(SIDE, G_FLAG, I_FLAG);

    ener = ising.energy(EXTFIELD);
    cout << "The energy is " << ener << endl;
    magn = ising.magnetization();
    cout << "The magnetization is " << magn << endl;
    ising.show_configuration();

    for(int i = 0; i < SUBSAMP; i++) ising.update(BETA, EXTFIELD);

    ener = ising.energy(EXTFIELD);
    cout << "The energy is " << ener << endl;
    magn = ising.magnetization();
    cout << "The magnetization is " << magn << endl;
    ising.show_configuration();

    ising.save_configuration();
}
