import os

import numpy as np
import matplotlib.pyplot as plt

os.system('g++ mersenne.cpp -o main')
os.system('./main')

a = np.loadtxt('random.dat', unpack='True')

#os.system('rm main random.dat')

x = a[1, 0:-2]
y = a[1, 1:-1]

plt.figure("Linear congruent generator")
plt.scatter(x, y, marker='+', label='Points')
plt.xlabel('$X_n$')
plt.ylabel('$X_{n+1}$')
plt.legend(loc='lower right')

x = a[2, 0:-2]
y = a[2, 1:-1]

plt.figure("Mersenne twister from C++")
plt.scatter(x, y, marker='+', label='Points')
plt.xlabel('$X_n$')
plt.ylabel('$X_{n+1}$')
plt.legend(loc='lower right')

plt.show()
