import os
import re

import numpy as np
import matplotlib.pyplot as plt



#-------------------------------------------------------------------------------

def plot_magnetization():
    plt.figure("fig_magnetization")
    plt.title("Average magnetization as a function of beta")
    plt.ylabel(r'$< M >$')
    plt.xlabel(r'$\beta$')

    sides = np.arange(20, 70, 10, dtype='int')

    for side in sides:
        x = []
        y = []
        directory = f"Side_{side}"
        print("logging: loading folder " + directory)
        for filename in os.listdir(directory):
            file = os.path.join(directory, filename)
            if os.path.isfile(file):
                data = np.loadtxt(file, unpack='True')
                beta = float(re.findall("\d+\.\d+", filename)[0])
                x.append(beta)
                y.append(abs(np.mean(data[1,:])))

        plt.scatter(x, y, s=10, label=f'side = {side}')

    plt.legend(loc='lower right')
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

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    # os.system('g++ ising_lattice.h ising_main.cpp -o main')
    # os.system('./main')
    # os.system('rm main')

    plot_magnetization()
