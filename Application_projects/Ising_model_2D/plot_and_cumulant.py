"""*****************************************************************************
*
* Plot program for the outcomes of the data analysis and fit Binder Cumulant
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
#*******************************************************************************

SIDE_MIN = 20
SIDE_MAX = 60
SIDE_SEP = 10

sides = np.arange(SIDE_MIN, SIDE_MAX+1, SIDE_SEP, dtype='int')

#--- Contents ------------------------------------------------------------------

def fit_fun(x, A, B):
    y = A * np.power(x, 2) + B
    return y

def load_data():
    """ Load data produced by analysis.cpp """

    data = {}
    for side in sides:
        # define data file path
        filename = f"side_{side}_data.dat"
        file_path = os.path.join("Data_analysis", filename)
        print("Loading " + file_path)
        # load data from each side file
        if os.path.isfile(file_path):
            data[side] = np.loadtxt(file_path, unpack='True')

    return data

def load_cumulant(data, n):
    """ Load cumulant data corresponding to betas[n] """

    recL = []
    cumulan = []
    cum_err = []
    for side in sides:
        # load cumulant data from analysis
        x, _, _, _, _, _, _, _, _, y, y_err = data[side]
        x, y, y_err = zip(*sorted(zip(x, y, y_err)))
        # organize data in function of side lenght
        recL.append(1 / side)
        cumulan.append(y[n])
        cum_err.append(y_err[n])

    return (x[n], recL, cumulan, cum_err)

def cumulant(data, n):
    """ Plot Binder cumulant as a function of L for n-th value of beta """

    #---Load points
    beta, recL, cumulan, cum_err = load_cumulant(data, n)

    #---Fit
    parameters, covariance = curve_fit(fit_fun, recL, cumulan, sigma=cum_err)
    fit_a = parameters[0]
    fit_b = parameters[1]
    std_deviation = np.sqrt(np.diag(covariance))
    fit_da = std_deviation[0]
    fit_db = std_deviation[1]
    print("\nFit parameters:")
    print(f"{fit_a} ± {fit_da}\n{fit_b} ± {fit_db}")
    # reduced chi squared
    fit_y = fit_fun(recL, fit_a, fit_b)
    chisq = np.sum(np.power(((cumulan - fit_y)/ cum_err), 2))
    chisqrd = chisq / (len(recL) - 3)
    print(f"Reduced chi squared: {chisqrd}")

    #---Plot
    beta = round(beta, 4)
    title = f"Binder cumulant beta = {beta}"
    print("\nPlot " + title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ C_B $')
    plt.xlabel(r'$ 1 / L $')
    # points and function
    fit_x = np.linspace(0., 0.051, 100)
    fit_y = fit_fun(fit_x, fit_a, fit_b)
    fit_label = f'fit {round(fit_b,4)} ± {round(fit_db,4)}'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(recL, cumulan, yerr=cum_err, fmt='<',label=f'beta: {beta}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_one(data):
    """ Plot one thing per time, e.g. susceptibility """

    title = "Plot susceptibility"
    print(title + "\n")
    # axis and style
    fig = plt.figure(title)
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \chi $')
    plt.xlabel(r'$\beta$')
    # load and plot susceptibility in function of beta
    for side in sides:
        x, _, _, _, _, _, _, y, y_err, _, _ = data[side]
        plt.errorbar(x, y, yerr=y_err, fmt='.', label=f'side = {side}')
    # save and show
    plt.legend(loc='upper right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_all(data):
    """ Plot everything all at once """

    title = "Plots from analysis"
    print(title + "\n")
    # axis and style
    fig, axes = plt.subplots(2, 2, num=title, figsize=(14, 14))
    plt.style.use('seaborn-whitegrid')
    # energy density
    axes[0, 0].set_title("Average energy density")
    axes[0, 0].set_ylabel(r'$ \langle \epsilon \rangle $')
    axes[0, 0].set_xlabel(r'$\beta$')
    # magnetization density
    axes[0, 1].set_title("Average magnetization density")
    axes[0, 1].set_ylabel(r'$\langle M \rangle $')
    axes[0, 1].set_xlabel(r'$\beta$')
    # specific heat
    axes[1, 0].set_title("Specific heat")
    axes[1, 0].set_ylabel(r'$ C_v$')
    axes[1, 0].set_xlabel(r'$\beta$')
    # susceptibility
    axes[1, 1].set_title("Susceptibility")
    axes[1, 1].set_ylabel(r'$ \chi $')
    axes[1, 1].set_xlabel(r'$\beta$')
    # load and plot data in function of beta
    for side in sides:
        betas, enes, enes_err, mags, mags_err, heat, heat_err, chi, chi_err, cum, cum_err = data[side]
        axes[0, 0].errorbar(betas, enes, yerr=enes_err, fmt='.', label=f'side = {side}')
        axes[0, 1].errorbar(betas, mags, yerr=mags_err, fmt='.', label=f'side = {side}')
        axes[1, 0].errorbar(betas, heat, yerr=heat_err, fmt='.', label=f'side = {side}')
        axes[1, 1].errorbar(betas, chi, yerr=chi_err, fmt='.', label=f'side = {side}')
    # legend, save and show
    axes[0, 0].legend(loc='lower left')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()

    plot_one(data)
    plot_all(data)

    #cumulant(data, 6)
    #cumulant(data, -10)
