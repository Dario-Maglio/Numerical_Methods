"""*****************************************************************************
*
* Plot and fit program for the outcomes of the data analysis
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt

#*******************************************************************************
# PARAMETERS OF THE SIMULATION
#
# SIDE_SEP = separation between the sides of different simulations.
#
# BETA_SEP = separation between the betas of different simulations.
#
#*******************************************************************************

SIDE_MIN = 20
SIDE_MAX = 60
SIDE_SEP = 10
BETA_INI = 0.3600
BETA_FIN = 0.5100
BETA_SEP = 0.0025

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')
betas = np.arange(BETA_INI, BETA_FIN, BETA_SEP, dtype='float')

#--- Contents ------------------------------------------------------------------

def load_data(data):
    """ Load data produced by analysis.cpp"""

    directory = "Data_analysis"

    for side in sides:
        filename = f"side_{side}_data.dat"
        file_path = os.path.join(directory, filename)
        print("Loading " + file_path)

        if os.path.isfile(file_path):
            data[side] = np.loadtxt(file_path, unpack='True')

def plot_one(data):
    """ Plot one thing per time """

    title = "Cumulante di Binder"
    fig = plt.figure(title)
    plt.title(title)
    plt.ylabel(r'$ C_B $')
    plt.xlabel(r'$\beta$')
    plt.style.use('seaborn-whitegrid')

    for side in sides:
        x, _, _, _, _, _, _, _, _, y, y_err = data[side]
        plt.errorbar(x, y, yerr=y_err, fmt='-x', label=f'side = {side}')

    print("\nPlot " + title + "\n")
    plt.legend(loc='lower left')
    plt.show()

def plot_all(data):
    """ Plot everything you can! """

    plt.style.use('seaborn-whitegrid')

    fig, axes = plt.subplots(2, 2, num="Plots from analysis", figsize=(14, 14))

    axes[0, 0].set_title("Average energy density")
    axes[0, 0].set_ylabel(r'$ \langle \epsilon \rangle $')
    axes[0, 0].set_xlabel(r'$\beta$')

    axes[0, 1].set_title("Average magnetization density")
    axes[0, 1].set_ylabel(r'$\langle M \rangle $')
    axes[0, 1].set_xlabel(r'$\beta$')

    axes[1, 0].set_title("Specific heat")
    axes[1, 0].set_ylabel(r'$ C_v$')
    axes[1, 0].set_xlabel(r'$\beta$')

    axes[1, 1].set_title("Susceptibility")
    axes[1, 1].set_ylabel(r'$ \chi $')
    axes[1, 1].set_xlabel(r'$\beta$')

    for side in sides:
        betas, enes, enes_err, mags, mags_err, heat, heat_err, chi, chi_err, cum, cum_err = data[side]

        axes[0, 0].errorbar(betas, enes, yerr=enes_err, fmt='.', label=f'side = {side}')

        axes[0, 1].errorbar(betas, mags, yerr=mags_err, fmt='.', label=f'side = {side}')

        axes[1, 0].errorbar(betas, heat, yerr=heat_err, fmt='.', label=f'side = {side}')

        axes[1, 1].errorbar(betas, chi, yerr=chi_err, fmt='.', label=f'side = {side}')


    print("\nPlots are ready! \n")
    axes[0, 0].legend(loc='lower left')
    plt.show()

def plot_cumulant(n, data):
    """ Plot Binder cumulant as a function of L for n-th value of beta """

    cum=[]
    cum_err=[]
    for side in sides:
        betas, _, _, _, _, _, _, _, _, y, y_err = data[side]
        cum.append(y[n])
        cum_err.append(y_err[n])

    title = f"Binder cumulant as a function of L for beta = {betas[n]}"
    fig = plt.figure(title)
    plt.title(title)

    plt.style.use('seaborn-whitegrid')
    plt.ylabel(r'$ C_B $')
    plt.xlabel(r'$ L $')

    plt.errorbar(sides, cum, yerr=cum_err, fmt='<',label=f'beta: {betas[n]}')

    print("\nPlot " + title + "\n")
    plt.legend(loc='lower left')
    plt.show()

def fit():
    pass

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    data = {}
    load_data(data)

    #plot_one(data)
    #plot_all(data)

    #plot_cumulant(3, data)
    plot_cumulant(52, data)

    #fit()
