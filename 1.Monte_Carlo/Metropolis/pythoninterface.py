import os

import numpy as np
import matplotlib.pyplot as plt

# os.system('g++ metropolis.cpp -o main')
# os.system('./main')
# os.system('rm main')

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------

y1 = y1[5000:15000]
print("10000 steps sample")
print(f"The average is {np.average(y1)} +- {np.sqrt(np.var(y1)/(len(y1) - 1))}")
print("Complete sample")
print(f"The average is {np.average(y2)} +- {np.sqrt(np.var(y2)/(len(y2) - 1))}")
