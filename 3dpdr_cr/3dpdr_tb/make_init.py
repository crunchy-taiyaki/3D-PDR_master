import numpy as np

mh = 1.67372E-24
RSUN = 6.957E10
MSUN = 1.988E33
G = 6.67259E-8

ncloud = 1E3
MCLOUD = 1E5
RCLOUD_cm = ((3.*MCLOUD*MSUN)/(4.*np.pi*ncloud*2.8*mh))**(1./3.)
RCLOUD_pc = RCLOUD_cm/(3.086E18)

#Set up the cloud
NRES = int(1E3)
R2 = np.log10(0.5)
arr1 = np.logspace(-3, R2, NRES/2, endpoint=False)
t1 = np.logspace(-3, R2, NRES/2, endpoint=True)
t1 -= 1.; t1 *= -1.
arr2 = t1[::-1]
rspace = np.concatenate((arr1, arr2))
rspace *= RCLOUD_pc

f = open("init_orion.dat", 'w')
for i in range(NRES):
	f.write("%e\t0.0\t0.0\t%e\n"%(rspace[i], 1E3))
f.close()
	
