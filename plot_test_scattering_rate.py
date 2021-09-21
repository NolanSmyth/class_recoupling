import numpy as np
import matplotlib.pyplot as plt

dat = np.genfromtxt("/Users/Nolan/Desktop/UCSC/Research/kineticRecoupling/class_kinetic_recoupling/data.csv", delimiter=",")
plt.plot(dat.T[0], dat.T[1])
plt.yscale('log')
plt.xscale('log')
plt.xlabel('T')
plt.ylabel('TSR')
plt.show()
