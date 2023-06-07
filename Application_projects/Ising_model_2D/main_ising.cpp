/*******************************************************************************
*
* Main program for the Ising simulations
*
*******************************************************************************/

//--- Preprocessor directives --------------------------------------------------

#include <iostream>
#include <chrono>

// Import the subroutine run_simulation
#include "ising_run_simulation.h"

using namespace std;


/*******************************************************************************
* PARAMETERS OF THE SIMULATION
*
* SIDE_SEP = separation between the sides of different simulations.
*
* BETA_SEP = separation between the betas of different simulations.
*
*******************************************************************************/

#define SIDE_SEP 10
#define SIDE_MIN 10
#define SIDE_MAX 70
// #define SIDE_MIN 70
// #define SIDE_MAX 71

// define outer betas -> 20 points
#define BETA_SEP 0.0050
#define BETA_INI 0.3600
#define BETA_FIN 0.4800
// define inner betas -> 50 points
#define BETA_C_SEP 0.00050
#define BETA_C_INI 0.41525
#define BETA_C_FIN 0.44000

//--- Main ---------------------------------------------------------------------

int main(){
    /* Main program iterates the timed Ising simulation over sides and betas */

    auto start = chrono::steady_clock::now();
    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
        for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
            run_simulation(side, beta);
        }
    }
    auto end = chrono::steady_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time: " << elapsed_seconds.count() << "s " << endl << endl;

    cout << "The work is done." << endl << endl;
    return 0;
}
