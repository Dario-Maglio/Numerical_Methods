"""*****************************************************************************
*
* Plot and fit program for the critical exponents
*
*****************************************************************************"""

import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat

#*******************************************************************************
# PARAMETERS OF THE SIMULATION
#
# SIDE_SEP = separation between the sides of different simulations.
#
#*******************************************************************************

SIDE_MIN = 20
SIDE_MAX = 60
SIDE_SEP = 10

sides = np.arange(SIDE_MIN, SIDE_MAX + 1, SIDE_SEP, dtype='int')

interval_chi = {20:(20,40), 30:(28,48), 40:(34,54), 50:(38,58), 60:(40,60)}
interval_cal = {20:(26,46), 30:(33,55), 40:(40,60), 50:(43,63), 60:(46,66)}
interval_mag = {20:(38,58), 30:(38,58), 40:(38,58), 50:(38,58), 60:(38,58)}

#--- Contents ------------------------------------------------------------------

def fit_par(x, a, b, c):
    y = a * np.power(x, 2) + b * np.power(x, 1) + c
    return y

def fit_beta(x, a, b, c):
    y = a / np.power(x, b) + c
    return y

def fit_chi(x, a, b):
    y = a * np.power(x, b)
    return y

def fit_mag(x, a, b):
    y = a / np.power(x, b)
    return y

def fit_cal(x, a, b):
    y = a / np.power(x, b)
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

#--- Parabolic fit -------------------------------------------------------------

def par_max_fit(fit_x, fit_y, fit_e):
    """ Get y_max values and correspinding x_max with parabolic fit """

    # fit the parabola
    parameters, covariance = curve_fit(fit_par, fit_x, fit_y, sigma=fit_e)
    std_deviation = np.sqrt(np.diag(covariance))
    # print values and uncertainties
    print("Fit parameters:")
    par_a = ufloat(parameters[0], std_deviation[0])
    print(par_a)
    par_b = ufloat(parameters[1], std_deviation[1])
    print(par_b)
    par_c = ufloat(parameters[2], std_deviation[2])
    print(par_c)
    # compute y_max and x_max
    print("Max point:")
    x_max = - par_b / (2 * par_a)
    print(x_max)
    y_max = (par_b * par_b)/(4 * par_a)
    y_max = y_max - (par_b * par_b)/(2 * par_a) + par_c
    print(y_max)
    # compute reduced chi squared
    fitted_y = fit_par(fit_x, *parameters)
    chisq = np.sum(np.power(((fit_y - fitted_y)/ fit_e), 2))
    chisqrd = chisq / (len(fit_x) - 4)
    print(f"Reduced chi squared: \n{chisqrd}\n")
    return (x_max, y_max, parameters)

def plot_par_chi(x, y, y_err, a, b, parameters, title):
    """ Plot parabolic fit for chi """

    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \chi $')
    plt.xlabel(r'$ \beta $')
    plt.xlim([x[a-5], x[b+5]])
    # plot data and fit in function of beta
    fit_x = np.linspace(x[a], x[b], 100)
    fit_y = fit_par(fit_x, *parameters)
    fit_label = f'fit [{round(x[a], 4)}, {round(x[b], 4)}]'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(x, y, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    path = os.path.join("Plots_and_fit", "Max_Sus")
    plt.savefig(os.path.join(path, title + ".png"))
    plt.show()

def plot_par_mag(x, y, y_err, a, b, parameters, title):
    """ Plot parabolic fit for mag """

    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle | M | \rangle $')
    plt.xlabel(r'$ \beta $')
    plt.xlim([x[a-5], x[b+5]])
    # plot data and fit in function of beta
    fit_x = np.linspace(x[a], x[b], 100)
    fit_y = fit_par(fit_x, *parameters)
    fit_label = f'fit [{round(x[a], 4)}, {round(x[b], 4)}]'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(x, y, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    path = os.path.join("Plots_and_fit", "Max_Mag")
    plt.savefig(os.path.join(path, title + ".png"))
    plt.show()

def plot_par_cal(x, y, y_err, a, b, parameters, title):
    """ Plot parabolic fit for cal """

    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ C_V $')
    plt.xlabel(r'$ \beta $')
    plt.xlim([x[a-5], x[b+5]])
    # plot data and fit in function of beta
    fit_x = np.linspace(x[a], x[b], 100)
    fit_y = fit_par(fit_x, *parameters)
    fit_label = f'fit [{round(x[a], 4)}, {round(x[b], 4)}]'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(x, y, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    path = os.path.join("Plots_and_fit", "Max_Cal")
    plt.savefig(os.path.join(path, title + ".png"))
    plt.show()

#--- Critical ratios and Beta --------------------------------------------------

def critical_point(beta_pc, beta_er):
    """ Get critical point and nu from beta pseudo critical """

    # exponential fit
    parameters, covariance = curve_fit(fit_beta, sides, beta_pc, sigma=beta_er)
    std_deviation = np.sqrt(np.diag(covariance))
    # print values and uncertainties
    print("Fit parameters:")
    par_a = ufloat(parameters[0], std_deviation[0])
    print(par_a)
    par_b = ufloat(parameters[1], std_deviation[1])
    print(par_b)
    par_c = ufloat(parameters[2], std_deviation[2])
    print(par_c)
    # compute reduced chi squared
    fitted_y = fit_beta(sides, *parameters)
    chisq = np.sum(np.power(((beta_pc - fitted_y)/ beta_er), 2))
    chisqrd = chisq / (len(sides) - 4)
    print(f"Reduced chi squared: \n{chisqrd}")
    # compute results
    beta_cr = par_c
    nu_exp = 1 / par_b
    return (beta_cr, nu_exp, parameters)

def critical_ratio(y_max, y_err, function):
    """ Get critical exponent ratio fitting a given function """

    # exponential fit
    parameters, covariance = curve_fit(function, sides, y_max, sigma=y_err)
    std_deviation = np.sqrt(np.diag(covariance))
    # print values and uncertainties
    print("Fit parameters:")
    par_a = ufloat(parameters[0], std_deviation[0])
    print(par_a)
    par_b = ufloat(parameters[1], std_deviation[1])
    print(par_b)
    # compute reduced chi squared
    fitted_y = function(sides, *parameters)
    chisq = np.sum(np.power(((y_max - fitted_y)/ y_err), 2))
    chisqrd = chisq / (len(sides) - 3)
    print(f"Reduced chi squared: \n{chisqrd}")
    # compute results
    ratio = par_b
    return (ratio, parameters)

#--- Plot subroutines ----------------------------------------------------------

def plot_beta_critical(beta_pc, beta_er, parameters):
    """ Plot beta_pc as a function of L """

    title = "Fit beta_pc as a function of L"
    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \beta_{pc} $')
    plt.xlabel(r'$ L $')
    # plot data and fit in function of beta
    fit_x = np.linspace(10, 80, 100)
    fit_y = fit_beta(fit_x, *parameters)
    fit_label = r'fit y = c + a / x^b'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(sides, beta_pc, yerr=beta_er, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_critical_chi(y_max, y_err, parameters):
    """ Plot chi_max as a function of L"""

    title = "Fit chi_max as a function of L"
    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \chi_{max} $')
    plt.xlabel(r'$ L $')
    # plot data and fit in function of beta
    fit_x = np.linspace(10, 80, 100)
    fit_y = fit_chi(fit_x, *parameters)
    fit_label = r'fit y = a * x^b'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(sides, y_max, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_critical_mag(y_max, y_err, parameters):
    """ Plot mag_max as a function of L"""

    title = "Fit mag_max as a function of L"
    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ \langle | M | \rangle $')
    plt.xlabel(r'$ L $')
    # plot data and fit in function of beta
    fit_x = np.linspace(10, 80, 100)
    fit_y = fit_mag(fit_x, *parameters)
    fit_label = r'fit y = a * x^b'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(sides, y_max, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

def plot_critical_cal(y_max, y_err, parameters):
    """ Plot cal_max as a function of L"""

    title = "Fit cal_max as a function of L"
    fig = plt.figure(title)
    # axis and style
    plt.style.use('seaborn-whitegrid')
    plt.title(title)
    plt.ylabel(r'$ C_V $')
    plt.xlabel(r'$ L $')
    # plot data and fit in function of beta
    fit_x = np.linspace(10, 80, 100)
    fit_y = fit_mag(fit_x, *parameters)
    fit_label = r'fit y = a * x^b'
    plt.plot(fit_x, fit_y, '-', label=fit_label)
    plt.errorbar(sides, y_max, yerr=y_err, fmt='.', label=f'simulation')
    # legend, save and show
    plt.legend(loc='lower right')
    plt.savefig(os.path.join("Plots_and_fit", title + ".png"))
    plt.show()

#--- Parabolic fit subroutines -------------------------------------------------

def parabolic_fit_chi(data):
    """ Subroutine for parabolic fit of chi max """

    beta_pc = []
    beta_er = []
    chi_max = []
    chi_err = []

    for side in sides:
        # select chi data to fit
        x, _, _, _, _, _, _, y, y_err, _, _ = data[side]
        x, y, y_err = zip(*sorted(zip(x, y, y_err)))
        a = interval_chi[side][0]
        b = interval_chi[side][1]
        # fit and append
        print(f"Susceptibility side {side}")
        x_max, y_max, parameters = par_max_fit(x[a:b], y[a:b], y_err[a:b])
        beta_pc.append(x_max.n)
        beta_er.append(x_max.std_dev)
        chi_max.append(y_max.n)
        chi_err.append(y_max.std_dev)
        # plot data and fit
        title = f"Fit max susceptibility side: {side} "
        print(title + "\n")
        plot_par_chi(x, y, y_err, a, b, parameters, title)

    return beta_pc, beta_er, chi_max, chi_err

def parabolic_fit_mag(data):
    """ Subroutine for parabolic fit of mag max """

    mag_max = []
    mag_err = []

    for side in sides:
        # select mag data to fit
        x, _, _, y, y_err, _, _, _, _, _, _ = data[side]
        x, y, y_err = zip(*sorted(zip(x, y, y_err)))
        a = interval_mag[side][0]
        b = interval_mag[side][1]
        # fit and append
        print(f"Magnetization side {side}")
        x_max, y_max, parameters = par_max_fit(x[a:b], y[a:b], y_err[a:b])
        mag_max.append(y_max.n)
        mag_err.append(y_max.std_dev)
        # plot data and fit
        title = f"Fit max magnetization side: {side} "
        print(title + "\n")
        plot_par_mag(x, y, y_err, a, b, parameters, title)

    return mag_max, mag_err

def parabolic_fit_cal(data):
    """ Subroutine for parabolic fit of cal max"""

    cal_max = []
    cal_err = []

    for side in sides:
        # select cal data to fit
        x, _, _, _, _, y, y_err, _, _, _, _ = data[side]
        x, y, y_err = zip(*sorted(zip(x, y, y_err)))
        a = interval_cal[side][0]
        b = interval_cal[side][1]
        # fit and append
        print(f"Specific heat side {side}")
        x_max, y_max, parameters = par_max_fit(x[a:b], y[a:b], y_err[a:b])
        cal_max.append(y_max.n)
        cal_err.append(y_max.std_dev)
        # plot data and fit
        title = f"Fit max specific heat side: {side} "
        print(title + "\n")
        plot_par_cal(x, y, y_err, a, b, parameters, title)

    return cal_max, cal_err

#--- Main ----------------------------------------------------------------------

if __name__ == '__main__':

    data = load_data()
    print("Loading complete! \n")

    #---Paraboli fits
    print("Start parabolic fit!\n")
    beta_pc, beta_er, chi_max, chi_err = parabolic_fit_chi(data)

    mag_max, mag_err = parabolic_fit_mag(data)

    cal_max, cal_err = parabolic_fit_cal(data)

    print("--- Study the critical point and exponents --------")

    #---Critical point
    print("Study pseudo critical beta")
    beta_cr, nu_exp, parameters = critical_point(beta_pc, beta_er)
    print("\nCritical beta:")
    print(beta_cr)
    print("Critical exponent nu:")
    print(nu_exp)

    print("\nPlot beta_pc as a function of L\n")
    plot_beta_critical(beta_pc, beta_er, parameters)

    #---Critical ratio gamma
    print("Study chi_max")
    ratio, parameters = critical_ratio(chi_max, chi_err, fit_chi)
    gamma_exp = nu_exp * ratio
    print("\nCritical exponent gamma:")
    print(gamma_exp)

    print("\nPlot chi_max as a function of L\n")
    plot_critical_chi(chi_max, chi_err, parameters)

    #---Critical ratio beta
    print("Study mag_max")
    ratio, parameters = critical_ratio(mag_max, mag_err, fit_mag)
    beta_exp = nu_exp * ratio
    print("\nCritical ratio:")
    print(ratio)
    print("Critical exponent beta:")
    print(beta_exp)

    print("\nPlot mag_max as a function of L\n")
    plot_critical_mag(mag_max, mag_err, parameters)

    #---Critical ratio alpha
    print("Study cal_max")
    ratio, parameters = critical_ratio(cal_max, cal_err, fit_cal)
    alpha_exp = nu_exp * ratio
    print("\nCritical exponent alpha:")
    print(alpha_exp)

    print("\nPlot chi_max as a function of L\n")
    plot_critical_cal(cal_max, cal_err, parameters)