from pylab import *
import numpy as np
import matplotlib.pyplot as plt


data=np.genfromtxt('B01.pdr.fin')
av_b01=data[:,2]
tgas_b01=data[:,3]
tdust_b01=data[:,4]
HI_b01=data[:,38]
H2_b01=data[:,39]
CII_b01=data[:,18]
CI_b01=data[:,32]
CO_b01=data[:,35]

data=np.genfromtxt('B02.pdr.fin')
av_b02=data[:,2]
tgas_b02=data[:,3]
tdust_b02=data[:,4]
HI_b02=data[:,38]
H2_b02=data[:,39]
CII_b02=data[:,18]
CI_b02=data[:,32]
CO_b02=data[:,35]

data=np.genfromtxt('B03.pdr.fin')
av_b03=data[:,2]
tgas_b03=data[:,3]
tdust_b03=data[:,4]
HI_b03=data[:,38]
H2_b03=data[:,39]
CII_b03=data[:,18]
CI_b03=data[:,32]
CO_b03=data[:,35]

data=np.genfromtxt('B04.pdr.fin')
av_b04=data[:,2]
tgas_b04=data[:,3]
tdust_b04=data[:,4]
HI_b04=data[:,38]
H2_b04=data[:,39]
CII_b04=data[:,18]
CI_b04=data[:,32]
CO_b04=data[:,35]

#Gas temperature
plt.figure()
plt.semilogx(av_b01,tgas_b01,label='B01',linewidth=2)
plt.semilogx(av_b02,tgas_b02,label='B02',linewidth=2)
plt.semilogx(av_b03,tgas_b03,label='B03',linewidth=2)
plt.semilogx(av_b04,tgas_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Gas temperature [K]',fontsize=15)
plt.grid()
savefig("Tgas.png",bbox_inches='tight')

#Dust temperature
plt.figure()
plt.semilogx(av_b01,tdust_b01,label='B01',linewidth=2)
plt.semilogx(av_b02,tdust_b02,label='B02',linewidth=2)
plt.semilogx(av_b03,tdust_b03,label='B03',linewidth=2)
plt.semilogx(av_b04,tdust_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Dust temperature [K]',fontsize=15)
plt.grid()
savefig("Tdust.png",bbox_inches='tight')

#Atomic hydrogen abundance
plt.figure()
plt.loglog(av_b01,HI_b01,label='B01',linewidth=2)
plt.loglog(av_b02,HI_b02,label='B02',linewidth=2)
plt.loglog(av_b03,HI_b03,label='B03',linewidth=2)
plt.loglog(av_b04,HI_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('HI Abundance',fontsize=15)
plt.grid()
savefig("HI.png",bbox_inches='tight')

#Molecular hydrogen abundance
plt.figure()
plt.loglog(av_b01,H2_b01,label='B01',linewidth=2)
plt.loglog(av_b02,H2_b02,label='B02',linewidth=2)
plt.loglog(av_b03,H2_b03,label='B03',linewidth=2)
plt.loglog(av_b04,H2_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'H$_2$ Abundance',fontsize=15)
plt.grid()
savefig("H2.png",bbox_inches='tight')

#Ionized carbon abundance
plt.figure()
plt.loglog(av_b01,CII_b01,label='B01',linewidth=2)
plt.loglog(av_b02,CII_b02,label='B02',linewidth=2)
plt.loglog(av_b03,CII_b03,label='B03',linewidth=2)
plt.loglog(av_b04,CII_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'CII Abundance',fontsize=15)
plt.grid()
savefig("CII.png",bbox_inches='tight')

#Atomic carbon abundance
plt.figure()
plt.loglog(av_b01,CI_b01,label='B01',linewidth=2)
plt.loglog(av_b02,CI_b02,label='B02',linewidth=2)
plt.loglog(av_b03,CI_b03,label='B03',linewidth=2)
plt.loglog(av_b04,CI_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'CI Abundance',fontsize=15)
plt.grid()
savefig("CI.png",bbox_inches='tight')

#Carbon monoxide abundance
plt.figure()
plt.loglog(av_b01,CO_b01,label='B01',linewidth=2)
plt.loglog(av_b02,CO_b02,label='B02',linewidth=2)
plt.loglog(av_b03,CO_b03,label='B03',linewidth=2)
plt.loglog(av_b04,CO_b04,label='B04',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'CO Abundance',fontsize=15)
plt.grid()
savefig("CO.png",bbox_inches='tight')


#plt.show()


