import os

import numpy as np
import matplotlib.pyplot as plt

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

#-------------------------------------------------------------------------------

def plot_ene_and_mag():
    sides = np.arange(SIDE_MIN, SIDE_MAX, SIDE_SEP, dtype='int')
    betas = np.arange(BETA_INI, BETA_FIN, BETA_SEP, dtype='float')

    fig, axes = plt.subplots(1, 2, num="energy and magnetization", figsize=(6, 6))

    axes[0].set_title("Average energy density")
    axes[0].set_ylabel(r'$< \epsilon >$')
    axes[0].set_xlabel(r'$\beta$')

    axes[1].set_title("Average magnetization density")
    axes[1].set_ylabel(r'$< M >$')
    axes[1].set_xlabel(r'$\beta$')

    for side in sides:
        ene = []
        mag = []

        directory = f"../Data_simulations/Side_{side}"
        print("Loading directory " + directory)

        for beta in betas:
            filename = "side_{0}_beta_{1:.6f}.dat".format(side, beta)
            file = os.path.join(directory, filename)

            if os.path.isfile(file):
                x, y = np.loadtxt(file, unpack='True')
                ene.append(np.mean(x))
                mag.append(np.mean(abs(y)))
        axes[0].plot(betas, ene, label=f'side = {side}')
        axes[1].plot(betas, mag, label=f'side = {side}')

    print("\nPlots of energy and magnetization: \n")
    axes[0].legend(loc='lower left')
    axes[1].legend(loc='lower right')
    plt.show()


def plot_all():
    plt.style.use('seaborn-whitegrid')

    fig, axes = plt.subplots(2, 2, num="specific heat and susceptibility", figsize=(14, 14))

    axes[0, 0].set_title("Average energy density")
    axes[0, 0].set_ylabel(r'$< \epsilon >$')
    axes[0, 0].set_xlabel(r'$\beta$')

    axes[0, 1].set_title("Average magnetization density")
    axes[0, 1].set_ylabel(r'$< M >$')
    axes[0, 1].set_xlabel(r'$\beta$')

    axes[1, 0].set_title("Specific heat")
    axes[1, 0].set_ylabel(r'$ V * \left( < \epsilon^2 > - < \epsilon >^2 \right)$')
    axes[1, 0].set_xlabel(r'$\beta$')

    axes[1, 1].set_title("Susceptibility")
    axes[1, 1].set_ylabel(r'$ V * \left( < M^2 > - < M >^2 \right)$')
    axes[1, 1].set_xlabel(r'$\beta$')

    sides = np.arange(20, 31, 10, dtype='int')

    directory = "../Data_analysis"


    for side in sides:
        filename = f"side_{side}_data.dat"
        file_path = os.path.join(directory, filename)
        print("Loading " + file_path)

        if os.path.isfile(file_path):
            betas, enes, enes_err, mags, mags_err, heat, heat_err, chi, chi_err, cum = np.loadtxt(file_path, unpack='True')

        axes[0, 0].errorbar(betas, enes, yerr=enes_err, fmt='.', label=f'side = {side}')

        axes[0, 1].errorbar(betas, mags, yerr=mags_err, fmt='.', label=f'side = {side}')

        axes[1, 0].errorbar(betas, heat, yerr=heat_err, fmt='.', label=f'side = {side}')

        axes[1, 1].errorbar(betas, chi, yerr=chi_err, fmt='.', label=f'side = {side}')


    print("\nPlots: \n")
    axes[0, 0].legend(loc='lower left')
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    # os.system('g++ ising_lattice.h ising_simulation.h main_ising.cpp -o main.out')
    # os.system('./main.out')
    # os.system('rm main.out')

    #plot_ene_and_mag()

    plot_all()
