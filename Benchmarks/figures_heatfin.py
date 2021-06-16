from pylab import *
import numpy as np
import matplotlib.pyplot as plt


data=np.genfromtxt('B01.heat.fin')
av_b01=data[:,2]
#hr01_b01=data[:,3]
hr02_b01=data[:,4]
#hr03_b01=data[:,5]
hr04_b01=data[:,6]
hr05_b01=data[:,7]
hr06_b01=data[:,8]
hr07_b01=data[:,9]
hr08_b01=data[:,10]
hr09_b01=data[:,11]
hr10_b01=data[:,12]
#hr11_b01=data[:,13]
hr12_b01=data[:,14]

data=np.genfromtxt('B02.heat.fin')
av_b02=data[:,2]
#hr01_b02=data[:,3]
hr02_b02=data[:,4]
#hr03_b02=data[:,5]
hr04_b02=data[:,6]
hr05_b02=data[:,7]
hr06_b02=data[:,8]
hr07_b02=data[:,9]
hr08_b02=data[:,10]
hr09_b02=data[:,11]
hr10_b02=data[:,12]
#hr11_b02=data[:,13]
hr12_b02=data[:,14]

data=np.genfromtxt('B03.heat.fin')
av_b03=data[:,2]
#hr01_b03=data[:,3]
hr02_b03=data[:,4]
#hr03_b03=data[:,5]
hr04_b03=data[:,6]
hr05_b03=data[:,7]
hr06_b03=data[:,8]
hr07_b03=data[:,9]
hr08_b03=data[:,10]
hr09_b03=data[:,11]
hr10_b03=data[:,12]
#hr11_b03=data[:,13]
hr12_b03=data[:,14]

data=np.genfromtxt('B04.heat.fin')
av_b04=data[:,2]
#hr01_b04=data[:,3]
hr02_b04=data[:,4]
#hr03_b04=data[:,5]
hr04_b04=data[:,6]
hr05_b04=data[:,7]
hr06_b04=data[:,8]
hr07_b04=data[:,9]
hr08_b04=data[:,10]
hr09_b04=data[:,11]
hr10_b04=data[:,12]
#hr11_b04=data[:,13]
hr12_b04=data[:,14]

#B01 
plt.figure()
plt.loglog(av_b01,hr02_b01,label='Photoelectric',linewidth=2)
plt.loglog(av_b01,hr04_b01,label='Carbon ionization',linewidth=2)
plt.loglog(av_b01,hr05_b01,label=r'H$_2$ formation',linewidth=2)
plt.loglog(av_b01,hr06_b01,label=r'H$_2$ photodiss.',linewidth=2)
plt.loglog(av_b01,hr07_b01,label='FUV pumbing',linewidth=2)
plt.loglog(av_b01,hr08_b01,label='Cosmic rays',linewidth=2)
plt.loglog(av_b01,hr09_b01,label='Microturbulent',linewidth=2)
plt.loglog(av_b01,hr10_b01,'--',label='Chemical',linewidth=2)
plt.loglog(av_b01,hr12_b01,'--',label='Gas-grain',linewidth=2)
plt.ylim([1e-42,1e-19])
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Heating function',fontsize=15)
plt.title('B01',fontsize=15,fontweight='bold')
plt.grid()
savefig("B01_heat.png",bbox_inches='tight')

#B02 
plt.figure()
plt.loglog(av_b02,hr02_b02,label='Photoelectric',linewidth=2)
plt.loglog(av_b02,hr04_b02,label='Carbon ionization',linewidth=2)
plt.loglog(av_b02,hr05_b02,label=r'H$_2$ formation',linewidth=2)
plt.loglog(av_b02,hr06_b02,label=r'H$_2$ photodiss.',linewidth=2)
plt.loglog(av_b02,hr07_b02,label='FUV pumbing',linewidth=2)
plt.loglog(av_b02,hr08_b02,label='Cosmic rays',linewidth=2)
plt.loglog(av_b02,hr09_b02,label='Microturbulent',linewidth=2)
plt.loglog(av_b02,hr10_b02,'--',label='Chemical',linewidth=2)
plt.loglog(av_b02,hr12_b02,'--',label='Gas-grain',linewidth=2)
plt.ylim([1e-42,1e-22])
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Heating function',fontsize=15)
plt.title('B02',fontsize=15,fontweight='bold')
plt.grid()
savefig("B02_heat.png",bbox_inches='tight')

#B03 
plt.figure()
plt.loglog(av_b03,hr02_b03,label='Photoelectric',linewidth=2)
plt.loglog(av_b03,hr04_b03,label='Carbon ionization',linewidth=2)
plt.loglog(av_b03,hr05_b03,label=r'H$_2$ formation',linewidth=2)
plt.loglog(av_b03,hr06_b03,label=r'H$_2$ photodiss.',linewidth=2)
plt.loglog(av_b03,hr07_b03,label='FUV pumbing',linewidth=2)
plt.loglog(av_b03,hr08_b03,label='Cosmic rays',linewidth=2)
plt.loglog(av_b03,hr09_b03,label='Microturbulent',linewidth=2)
plt.loglog(av_b03,hr10_b03,'--',label='Chemical',linewidth=2)
plt.loglog(av_b03,hr12_b03,'--',label='Gas-grain',linewidth=2)
plt.ylim([1e-34,1e-19])
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Heating function',fontsize=15)
plt.title('B03',fontsize=15,fontweight='bold')
plt.grid()
savefig("B03_heat.png",bbox_inches='tight')

#B04 
plt.figure()
plt.loglog(av_b04,hr02_b04,label='Photoelectric',linewidth=2)
plt.loglog(av_b04,hr04_b04,label='Carbon ionization',linewidth=2)
plt.loglog(av_b04,hr05_b04,label=r'H$_2$ formation',linewidth=2)
plt.loglog(av_b04,hr06_b04,label=r'H$_2$ photodiss.',linewidth=2)
plt.loglog(av_b04,hr07_b04,label='FUV pumbing',linewidth=2)
plt.loglog(av_b04,hr08_b04,label='Cosmic rays',linewidth=2)
plt.loglog(av_b04,hr09_b04,label='Microturbulent',linewidth=2)
plt.loglog(av_b04,hr10_b04,'--',label='Chemical',linewidth=2)
plt.loglog(av_b04,hr12_b04,'--',label='Gas-grain',linewidth=2)
plt.ylim([1e-32,1e-18])
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Heating function',fontsize=15)
plt.title('B04',fontsize=15,fontweight='bold')
plt.grid()
savefig("B04_heat.png",bbox_inches='tight')

#plt.show()
