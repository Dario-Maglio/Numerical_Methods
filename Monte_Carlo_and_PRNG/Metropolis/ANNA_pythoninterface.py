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

def plot_box_muller():
    a = np.loadtxt('file_box_muller.dat', unpack='True' )
    x = a[0, :]
    y = a[1, :]

    plt.figure('Box Muller Guassian Distribution')
    plt.scatter(x, y, label='Points')
    plt.xlabel('Steps')
    plt.ylabel('Variable')

    plt.figure("BM histogram gaussian")
    plt.hist(y, label='BM points', bins=300)
    plt.ylabel('Probability density function')
    plt.xlabel('Domain')
    plt.legend(loc='upper right')
    
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

    fig, axes = plt.subplots(1, 2, sharex=True, num="MC gaussian autoc")

    axes[0].set_title("Autocorrelation \n k = |j - i|")
    axes[0].plot(x, c, marker='.', markersize=0.1)
    axes[0].set_ylabel('C(k)')

    axes[1].set_title("Integrated autocorrelation time \n dim sample = 2^22")
    axes[1].plot(x, t, marker='.', markersize=0.1)
    axes[1].set_ylabel('tau_int(k)')

    plt.show()

def plot_bootstrap():
    a = np.loadtxt('file_boot_w.dat', unpack='True')
    xw = a[0, :]
    yw = a[1, :]
    a = np.loadtxt('file_boot_c.dat', unpack='True')
    xc = a[0, :]
    yc = a[1, :]

    txt = (" samp 300 -> 0.997062 ± 0.00159631 |"
        + " lenght 1 -> 0.99708 ± 0.00159522 \n"
        + " lenght 65536 -> 0.997017 ± 0.00257424 |"
        + " lenght 2048 -> 0.997222 ± 0.0028662")


    fig = plt.figure("MC bootstrap")
    fig.suptitle("Bootstrap for the estimator \n" +
                 r"$\dfrac{< x^4 >}{3 < x^2 >^2}$" +
                 r"   with dim samples = $10^6$")

    ax = fig.add_subplot(121)
    ax.set_title("without correlations \n lenght correlated blocks = 1")
    ax.scatter(xw, yw, marker='+')
    ax.set_ylabel('sigma')
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.set_xlabel('number fake samples')


    ax = fig.add_subplot(122)
    ax.set_title("with correlations \n num fake samples = 300")
    ax.scatter(xc, yc, marker='+')
    ax.set_ylabel('sigma')
    ax.set_xlabel('lenght correlated blocks')
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.set_xscale('log', base=2)

    fig.text(.5, .03, txt, ha='center')

    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    # os.system('g++ metropolis.cpp -o main')
    # os.system('./main')
    # os.system('rm main')

    # plot_metropolis()

    # plot_averages()

    # plot_autocorr()

    # plot_bootstrap()
    
    plot_box_muller()


