import matplotlib.pyplot as plt
import numpy as np
import sys

dat = np.loadtxt("TAURUS2.pdr.fin")
zeta_dat = np.loadtxt("zeta.txt")

av = dat[:,2]
temp = dat[:,3]
crir = zeta_dat[:,1]*1.3E-17

print len(av)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.loglog(av, temp, 'k-', linewidth=2,label="Temp")
ax2.loglog(av, crir, 'r--', linewidth=2, label="CRIR")
ax.set_xlabel("A$_V$", fontsize=16)
ax.set_ylabel("Temperature", fontsize=16)
ax2.set_ylabel("CRIR", fontsize=16, color='red')
ax.set_xlim(min(av), max(av))
fig.savefig("temperature.png")
plt.show()
