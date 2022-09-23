import os

import numpy as np
import matplotlib.pyplot as plt

os.system('g++ metropolis.cpp -o main')
os.system('./main')

a = np.loadtxt('f_metropolis.dat', unpack='True')

os.system('rm main')

x = a[0, :]
y = a[1, :]

plt.figure("MC evolution gaussian")
plt.plot(x, y, marker='+', markersize=0.1, label='MC points')
plt.xlabel('Evolution')
plt.ylabel('Distribution')
plt.legend(loc='lower right')

y = y[15000:]

plt.figure("MC histogram gaussian")
plt.hist(y, label='MC points', bins=500)
plt.ylabel('Probability density function')
plt.xlabel('Domain')
plt.legend(loc='upper right')

plt.show()
