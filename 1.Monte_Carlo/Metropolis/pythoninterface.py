import os

import numpy as np
import matplotlib.pyplot as plt

THERMALIZATION = 15000

#-------------------------------------------------------------------------------

def plot_metropolis():
    a = np.loadtxt('file_metropol.dat', unpack='True')

    x = [i for i in range(THERMALIZATION)]
    y1 = a[:THERMALIZATION]
    y2 = a[THERMALIZATION:]

    fig, axes = plt.subplots(1, 2, num="MC gaussian", figsize=(12, 12))

    axes[0].set_title("Thermalization phase \n sequence of 15000 steps")
    axes[0].plot(x, y1, marker='.', markersize=0.01,
                 label='start = 0 , delta = 0.1')
    axes[0].set_ylabel('Evolution')
    axes[0].set_xlabel('Step')
    axes[0].legend(loc='lower right')

    axes[1].set_title("Histogram after thermalization \n aver = 5 , sigma = 1")
    axes[1].hist(y2, bins=500, density=True)
    axes[1].set_ylabel('Counting')
    axes[1].set_xlabel('Variable')

    plt.show()

def plot_averages():
    a = np.loadtxt('file_averages.dat', unpack='True')
    x = [i for i in range(len(a))]

    plt.figure("MC gaussian averages")
    plt.title("Sample averages \n dim subsamples = 10000")
    plt.scatter(x, a, marker='.')
    plt.ylim([0, 7])
    plt.ylabel('Average')
    plt.xlabel('Sample')

    plt.show()

def plot_autocorr():
    a = np.loadtxt('file_autocor.dat', unpack='True')
    x = a[0, :]
    c = a[1, :]
    t = a[2, :]

    fig, axes = plt.subplots(1, 2, sharex=True, num="MC autocorrelation")

    axes[0].set_title("Autocorrelation \n k = |j - i|")
    axes[0].plot(x, c, marker='.', markersize=0.1)
    axes[0].set_ylabel('C(k)')

    axes[1].set_title("Integrated autocorrelation time \n dim sample = 2^22")
    axes[1].plot(x, t, marker='.', markersize=0.1)
    axes[1].set_ylabel('tau_int(k)')

    plt.show()

def plot_bootstrap():
    a = np.loadtxt('file_boot.dat', unpack='True')
    x = a[0, :]
    y = a[1, :]

    plt.figure("MC bootstrap without corr")
    plt.title("Bootstrap \n k = number of fake samples, dim samples = 2^22")
    plt.plot(x, y, marker='.', markersize=0.1)
    plt.ylabel('sigma(k)')
    plt.xlabel('step k')

    plt.show()

def plot_bootstrap_corr():
    a = np.loadtxt('file_boot_corr.dat', unpack='True')
    x = a[0, :]
    y = a[1, :]

    plt.figure("MC bootstrap with corr")
    plt.title("Bootstrap with correlations \n" +
              "correl block dim = 2^k, dim sample = 2^22, num samples = 100")
    plt.plot(x, y, marker='.', markersize=0.1)
    plt.ylabel('sigma(k)')
    plt.xlabel('step k')

    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    # os.system('g++ metropolis.cpp -o main')
    # os.system('./main')
    # os.system('rm main')

    plot_metropolis()

    plot_averages()

    plot_autocorr()

    plot_bootstrap()

    plot_bootstrap_corr()
