"""*****************************************************************************
*
* Plot and fit program for the outcomes of the data analysis
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

    title = "Susceptibility"
    print("\nPlot " + title + "\n")

    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \chi $')
    plt.xlabel(r'$\beta$')

    for side in sides:
        """ betas, enes, enes_err, mags, mags_err,
        heat, heat_err, chi, chi_err, cum, cum_err """
        x, _, _, _, _, _, _, y, y_err, _, _ = data[side]
        plt.errorbar(x, y, yerr=y_err, fmt='.', label=f'side = {side}')

    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_all(data):
    """ Plot everything you can! """

    print("\nPlot simulation! \n")
    fig, axes = plt.subplots(2, 2, num="Plots from analysis", figsize=(14, 14))
    plt.style.use('seaborn-whitegrid')

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

    axes[0, 0].legend(loc='lower left')
    plt.savefig(os.path.join("Plots_and_fit", "Simulation.png"))
    plt.show()

def function(x, a, b):
    y = a*np.power(x, 2) + b
    return y

def cumulant(n, data):
    """ Plot Binder cumulant as a function of L for n-th value of beta """

    x = []
    cumulan = []
    cum_err = []
    for side in sides:
        betas, _, _, _, _, _, _, _, _, y, y_err = data[side]
        x.append(1 / side)
        cumulan.append(y[n])
        cum_err.append(y_err[n])

    # Fit
    parameters, covariance = curve_fit(function, x, cumulan, sigma=cum_err)
    fit_a = parameters[0]
    fit_b = parameters[1]
    std_deviation = np.sqrt(np.diag(covariance))
    fit_da = std_deviation[0]
    fit_db = std_deviation[1]
    print("\nFit parameters:")
    print(f"{fit_a} \pm {fit_da} | {fit_b} \pm {fit_db}")

    chisq = np.sum(np.power((cumulan - function(x, fit_a, fit_b)/ cum_err), 2))
    chisqrd = chisq / (len(x)-3)
    print(f"Reduced chi squared: {chisq}")

    # Plot
    title = f"Binder cumulant beta = {betas[n]}"
    print("\nPlot " + title + "\n")

    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ C_B $')
    plt.xlabel(r'$ 1 / L $')

    fit_x = np.linspace(0., 0.051, 100)
    plt.plot(fit_x, function(fit_x, fit_a, fit_b), '-', label='fit')
    plt.errorbar(x, cumulan, yerr=cum_err, fmt='<',label=f'beta: {betas[n]}')

    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    data = {}
    load_data(data)

    #plot_one(data)
    #plot_all(data)

    cumulant(8, data)
    cumulant(50, data)

    #popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
