import os

import numpy as np
import matplotlib.pyplot as plt

os.system('g++ metropolis.cpp -o main')
os.system('./main')

a = np.loadtxt('f_metropolis.dat', unpack='True')

#os.system('rm main random.dat')

x = a[0, 10000:]
y = a[1, 10000:]

plt.figure("MC evolution gaussian")
plt.plot(x, y, marker='+', label='MC points')
plt.xlabel('Evolution')
plt.ylabel('Distribution')
plt.legend(loc='lower right')

plt.figure("MC histogram gaussian")
plt.hist(y, label='MC points', bins=300)
plt.ylabel('Probability density function')
plt.xlabel('Domain')
plt.legend(loc='upper right')

plt.show()
