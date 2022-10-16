import os

import numpy as np
import matplotlib.pyplot as plt



#-------------------------------------------------------------------------------

def plot_magnetization():
    plt.figure("magnetization")
    plt.title("Average magnetization as a function of beta")
    plt.ylabel(r'$< M >$')
    plt.xlabel(r'$\beta$')

    sides = np.arange(20, 61, 10, dtype='int')
    betas = np.arange(0.3600, 0.5101, 0.0025, dtype='float')

    for side in sides:
        ene = []
        mag = []

        directory = f"Side_{side}"
        print("logging: loading directory " + directory)

        for beta in betas:
            filename = "side_{0}_beta_{1:.6f}.dat".format(side, beta)
            file = os.path.join(directory, filename)
            print("logging: loading file -> " + file)
            if os.path.isfile(file):
                x, y = np.loadtxt(file, unpack='True')
                ene.append(np.mean(x))
                mag.append(np.mean(abs(y)))

        plt.plot(betas, mag, label=f'side = {side}')

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
