/*******************************************************************************
*
* Main program for the parallel Ising simulations
*
*******************************************************************************/

//----Preprocessor directives---------------------------------------------------
#include <iostream>
#include <vector>
#include <thread>
#include <chrono>

// Import the simulation subroutine
#include "ising_run_simulation.h"

/*
* CONFIGURATION PARAMETERS OF THE SIMULATION
* BETA_SEP = separation between the betas of different simulations.
* SIDE_SEP = separation between the sides of different simulations.
*/
#define BETA_INI 0.3600
#define BETA_FIN 0.5100
#define BETA_SEP 0.0025
#define SIDE_MIN 20
#define SIDE_MAX 60
#define SIDE_SEP 10

using namespace std;

//----Contents------------------------------------------------------------------


/* Main program iterates over sides and betas */
int main(){
    vector<thread> threadPool;

    auto start = chrono::steady_clock::now();
    for(int side = SIDE_MIN; side <= SIDE_MAX; side += SIDE_SEP){
      for(float beta = BETA_INI; beta <= BETA_FIN; beta += BETA_SEP){
        threadPool.emplace_back([side, beta]() {run_simulation(side, beta); });
      }
    }

    for (auto& thr : threadPool) thr.join();
    auto end = chrono::steady_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Elapsed time: " << elapsed_seconds.count() << "s " << endl << endl;

    cout << "The work is done." << endl << endl;
    return 0;
}
