import os

import numpy as np
import matplotlib.pyplot as plt

os.system('gfortran library.f90 Test_drive_basic.f90 -o drive.x')
os.system('./drive.x')

a = np.loadtxt('function.dat', unpack='True')

os.system('rm drive.x function.dat library.mod')

x = a[0, :]
y = a[2, :]

plt.figure("Figura 1")
plt.errorbar(x, y, label='Points')
plt.xlabel('Numeri')
plt.ylabel('Risultati')
plt.legend(loc='lower right')
plt.show()
