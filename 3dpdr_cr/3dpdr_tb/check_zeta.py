import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt("zeta_test.txt")
plt.figure()
plt.semilogy(dat[:,0], dat[:,1]*1.3E-17)
plt.xlabel(r"X (pc)", fontsize=16)
plt.ylabel(r"$\zeta$ (s$^{-1}$)", fontsize=16)
plt.xlim(dat[:,0][0], dat[:,0][-1])
plt.show()

R = 3.228
dat2 = np.loadtxt("colTest.dat")
xm = dat2[:,1]
x = dat2[:,0]

plt.plot(x, xm)
plt.show()
