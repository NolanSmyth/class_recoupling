import matplotlib.pyplot as plt
import numpy as np
import glob

datdir = "./data/"

filenames = glob.glob(datdir + "phase*_pk.dat")

filename = filenames[0]

data = np.loadtxt(filename[0])
plt.plot(data[:, 0], data[:, 1])
plt.yscale("log")
plt.xscale("log")
plt.xlim(1, 1e2)
# plt.savefig(filename + ".png")
plt.show()
