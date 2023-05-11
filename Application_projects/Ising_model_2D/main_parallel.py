"""/****************************************************************************
*
* Main program for the parallel Ising simulations
*
****************************************************************************/"""

#-----Initial imports and definitions-------------------------------------------

import time

import ROOT
import numpy as np
import multiprocessing as mp

ROOT.gInterpreter.ProcessLine('#include "ising_lattice_class.h"')
ROOT.gInterpreter.ProcessLine('#include "ising_run_simulation.h"')

"""
* CONFIGURATION PARAMETERS OF THE LATTICE
* BETA_SEP = separation between the betas of different simulations.
* SIDE_SEP = separation between the sides of different simulations.
"""
BETA_INI = 0.3600
BETA_FIN = 0.5101
BETA_SEP = 0.0025
SIDE_MIN = 20
SIDE_MAX = 61
SIDE_SEP = 10

sides = np.arange(SIDE_MIN, SIDE_MAX, SIDE_SEP, dtype='int')
betas = np.arange(BETA_INI, BETA_FIN, BETA_SEP, dtype='float')
args = np.meshgrid(sides, betas)
sides, betas = np.reshape(args, (2, -1))

#-----Contents------------------------------------------------------------------

if __name__ == '__main__':
    # Define an output queue if you need results
    output = mp.Queue()

    # Setup a list of processes that we want to run
    processes = [
        mp.Process(
            target=ROOT.run_simulation,
            args=(int(side), beta)
        ) for side, beta in zip(sides, betas)
    ]

    # Clocking processes
    tic = time.perf_counter()
    # Run processes
    for p in processes:
        p.start()
    # Exit the completed processes
    for p in processes:
        p.join()
    toc = time.perf_counter()

    print(f"Executed in {toc - tic:0.4f} seconds")

    print("The work is done.")
