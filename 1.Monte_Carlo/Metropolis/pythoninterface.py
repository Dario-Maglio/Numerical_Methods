import os

import numpy as np
import matplotlib.pyplot as plt



#-------------------------------------------------------------------------------

def plot_metropolis():
    a = np.loadtxt('f_metropolis.dat', unpack='True')

    x = [i+1 for i in range(15000)]
    y1 = a[:15000]
    y2 = a[15000:]

    plt.figure("MC evolution gaussian")
    plt.title("Thermalization" + "\n" + "Sequence of 15000 Metropolis steps")
    plt.plot(x, y1, marker='.', markersize=0.1, label='MC points')
    plt.xlabel('Step')
    plt.ylabel('Evolution')
    plt.legend(loc='lower right')

    plt.figure("MC histogram gaussian")
    plt.title("Histogram of the MC steps" + "\n" + "average = 5. , sigma = 1.")
    plt.hist(y2, bins=500)
    plt.ylabel('Counting')
    plt.xlabel('MC step')

    plt.show()

def plot_averages():
    a = np.loadtxt('f_averages.dat', unpack='True')
    x = [i+1 for i in range(len(a))]

    a = a[15:]
    x = x[15:]

    plt.figure("MC averages gaussian")
    plt.title("Scatter of the sample averages")
    plt.scatter(x, a, marker='.')
    plt.ylabel('Average')
    plt.ylim([0, 7.1])
    plt.xlabel('Sample')

    plt.show()

def plot_autocorr():
    a = np.loadtxt('f_autocorr.dat', unpack='True')
    x = [i+1 for i in range(len(a))]

    plt.figure("MC autocorrelation")
    plt.title("Autocorrelation of a 1000000 steps sample")
    plt.plot(x, a, marker='.', markersize=0.1)
    plt.ylabel('C(k)')
    plt.xlabel('k')

    a = np.loadtxt('f_tau_int.dat', unpack='True')
    x = [i+1 for i in range(len(a))]

    plt.figure("MC autocorrelation time")
    plt.title("Integrated autocorrelation time")
    plt.plot(x, a, marker='.', markersize=0.1)
    plt.ylabel('tau_int(k)')
    plt.xlabel('k')

    plt.show()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    os.system('g++ metropolis.cpp -o main')
    os.system('./main')
    os.system('rm main')

    # plot_metropolis()

    # plot_averages()

    plot_autocorr()
