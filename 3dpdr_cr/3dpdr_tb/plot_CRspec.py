import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl

label_size=16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

xarr = []

f = open(sys.argv[1], 'r')
fileI = f.readlines()
f.close()
dat = fileI[0].split()
Nene = int(dat[0])
Earr = np.zeros(Nene, 'float64')
Specs = []
xs = []
i = 0
for d in dat[1:]:
	Earr[i] = float(d)
	i += 1

i = 0
for l in fileI[1:]:
	dat = l.split()
	xs.append(float(dat[1]))
	speci = np.zeros(Nene, 'float64')
	j = 0
	for d in dat[2:]:
		speci[j] = float(d)
		j += 1
	Specs.append(speci)
	i += 1

minX = min(xs); maxX = max(xs)
cmap = plt.get_cmap('viridis')
cNorm = colors.Normalize(vmin=np.log10(minX), vmax=np.log10(maxX))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
scalarMap._A = []

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
for i in range(len(xs)):
	lx = np.log10(xs[i])
	color = scalarMap.to_rgba(lx)
	ax.loglog(Earr, Specs[i], linewidth=2, color=color)
ax.set_xlabel("E (eV)", fontsize=16)
ax.set_ylabel(r"$j_p(E)$ (protons cm$^{-2}$ s$^{-1}$ eV$^{-1}$ sr$^{-1}$)", fontsize=16)
c1 = fig.colorbar(scalarMap, ax=ax)
c1.set_label(r"$\log x$ (pc)", fontsize=14)
plt.tight_layout()
fig.savefig("Example_spec.png")
plt.show()
