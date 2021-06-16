from pylab import *
import numpy as np
import matplotlib.pyplot as plt


data=np.genfromtxt('B01.line.fin')
av_b01=data[:,2]
l01_b01=data[:,3] #CII
l02_b01=data[:,4] #CI 1-0
l03_b01=data[:,5] #CI 2-0
l04_b01=data[:,6] #CI 2-1
l05_b01=data[:,7] #OI 1-0
l06_b01=data[:,8] #OI 2-0
l07_b01=data[:,9] #OI 2-1
l08_b01=data[:,10] #CO (1-0)
l09_b01=data[:,11] #CO (2-1)
l10_b01=data[:,12] #CO (3-2)
l11_b01=data[:,13] #CO (4-3)
l12_b01=data[:,14] #CO (5-4)
l13_b01=data[:,15] #CO (6-5)
l14_b01=data[:,16] #CO (7-6)
l15_b01=data[:,17] #CO (8-7)
l16_b01=data[:,18] #CO (9-8)
l17_b01=data[:,19] #CO (10-9)

data=np.genfromtxt('B02.line.fin')
av_b02=data[:,2]
l01_b02=data[:,3] #CII
l02_b02=data[:,4] #CI 1-0
l03_b02=data[:,5] #CI 2-0
l04_b02=data[:,6] #CI 2-1
l05_b02=data[:,7] #OI 1-0
l06_b02=data[:,8] #OI 2-0
l07_b02=data[:,9] #OI 2-1
l08_b02=data[:,10] #CO (1-0)
l09_b02=data[:,11] #CO (2-1)
l10_b02=data[:,12] #CO (3-2)
l11_b02=data[:,13] #CO (4-3)
l12_b02=data[:,14] #CO (5-4)
l13_b02=data[:,15] #CO (6-5)
l14_b02=data[:,16] #CO (7-6)
l15_b02=data[:,17] #CO (8-7)
l16_b02=data[:,18] #CO (9-8)
l17_b02=data[:,19] #CO (10-9)

data=np.genfromtxt('B03.line.fin')
av_b03=data[:,2]
l01_b03=data[:,3] #CII
l02_b03=data[:,4] #CI 1-0
l03_b03=data[:,5] #CI 2-0
l04_b03=data[:,6] #CI 2-1
l05_b03=data[:,7] #OI 1-0
l06_b03=data[:,8] #OI 2-0
l07_b03=data[:,9] #OI 2-1
l08_b03=data[:,10] #CO (1-0)
l09_b03=data[:,11] #CO (2-1)
l10_b03=data[:,12] #CO (3-2)
l11_b03=data[:,13] #CO (4-3)
l12_b03=data[:,14] #CO (5-4)
l13_b03=data[:,15] #CO (6-5)
l14_b03=data[:,16] #CO (7-6)
l15_b03=data[:,17] #CO (8-7)
l16_b03=data[:,18] #CO (9-8)
l17_b03=data[:,19] #CO (10-9)

data=np.genfromtxt('B04.line.fin')
av_b04=data[:,2]
l01_b04=data[:,3] #CII
l02_b04=data[:,4] #CI 1-0
l03_b04=data[:,5] #CI 2-0
l04_b04=data[:,6] #CI 2-1
l05_b04=data[:,7] #OI 1-0
l06_b04=data[:,8] #OI 2-0
l07_b04=data[:,9] #OI 2-1
l08_b04=data[:,10] #CO (1-0)
l09_b04=data[:,11] #CO (2-1)
l10_b04=data[:,12] #CO (3-2)
l11_b04=data[:,13] #CO (4-3)
l12_b04=data[:,14] #CO (5-4)
l13_b04=data[:,15] #CO (6-5)
l14_b04=data[:,16] #CO (7-6)
l15_b04=data[:,17] #CO (8-7)
l16_b04=data[:,18] #CO (9-8)
l17_b04=data[:,19] #CO (10-9)


#B01-fine structure
plt.figure()
plt.loglog(av_b01,l01_b01,label=r'CII 158$\mu$m',linewidth=2)
plt.loglog(av_b01,l02_b01,label=r'CI 609$\mu$m',linewidth=2)
plt.loglog(av_b01,l04_b01,label=r'CI 320$\mu$m',linewidth=2)
plt.loglog(av_b01,l05_b01,label=r'OI 64$\mu$m',linewidth=2)
plt.loglog(av_b01,l07_b01,label=r'OI 146$\mu$m',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B01 Fine structure',fontsize=15,fontweight='bold')
plt.grid()
savefig("B01_fine.png",bbox_inches='tight')
#B01-CO ladder
plt.figure()
plt.loglog(av_b01,l08_b01,label=r'CO(1-0)',linewidth=2)
plt.loglog(av_b01,l09_b01,label=r'CO(2-1)',linewidth=2)
plt.loglog(av_b01,l10_b01,label=r'CO(3-2)',linewidth=2)
plt.loglog(av_b01,l11_b01,label=r'CO(4-3)',linewidth=2)
plt.loglog(av_b01,l12_b01,label=r'CO(5-4)',linewidth=2)
plt.loglog(av_b01,l13_b01,label=r'CO(6-5)',linewidth=2)
plt.loglog(av_b01,l14_b01,label=r'CO(7-6)',linewidth=2)
plt.loglog(av_b01,l15_b01,'--',label=r'CO(8-7)',linewidth=2)
plt.loglog(av_b01,l16_b01,'--',label=r'CO(9-8)',linewidth=2)
plt.loglog(av_b01,l17_b01,'--',label=r'CO(10-9)',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B01 CO ladder',fontsize=15,fontweight='bold')
plt.grid()
savefig("B01_ladder.png",bbox_inches='tight')

#B02-fine structure
plt.figure()
plt.loglog(av_b02,l01_b02,label=r'CII 158$\mu$m',linewidth=2)
plt.loglog(av_b02,l02_b02,label=r'CI 609$\mu$m',linewidth=2)
plt.loglog(av_b02,l04_b02,label=r'CI 320$\mu$m',linewidth=2)
plt.loglog(av_b02,l05_b02,label=r'OI 64$\mu$m',linewidth=2)
plt.loglog(av_b02,l07_b02,label=r'OI 146$\mu$m',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B02 Fine structure',fontsize=15,fontweight='bold')
plt.grid()
savefig("B02_fine.png",bbox_inches='tight')
#B02-CO ladder
plt.figure()
plt.loglog(av_b02,l08_b02,label=r'CO(1-0)',linewidth=2)
plt.loglog(av_b02,l09_b02,label=r'CO(2-1)',linewidth=2)
plt.loglog(av_b02,l10_b02,label=r'CO(3-2)',linewidth=2)
plt.loglog(av_b02,l11_b02,label=r'CO(4-3)',linewidth=2)
plt.loglog(av_b02,l12_b02,label=r'CO(5-4)',linewidth=2)
plt.loglog(av_b02,l13_b02,label=r'CO(6-5)',linewidth=2)
plt.loglog(av_b02,l14_b02,label=r'CO(7-6)',linewidth=2)
plt.loglog(av_b02,l15_b02,'--',label=r'CO(8-7)',linewidth=2)
plt.loglog(av_b02,l16_b02,'--',label=r'CO(9-8)',linewidth=2)
plt.loglog(av_b02,l17_b02,'--',label=r'CO(10-9)',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B02 CO ladder',fontsize=15,fontweight='bold')
plt.grid()
savefig("B02_ladder.png",bbox_inches='tight')

#B03-fine structure
plt.figure()
plt.loglog(av_b03,l01_b03,label=r'CII 158$\mu$m',linewidth=2)
plt.loglog(av_b03,l02_b03,label=r'CI 609$\mu$m',linewidth=2)
plt.loglog(av_b03,l04_b03,label=r'CI 320$\mu$m',linewidth=2)
plt.loglog(av_b03,l05_b03,label=r'OI 64$\mu$m',linewidth=2)
plt.loglog(av_b03,l07_b03,label=r'OI 146$\mu$m',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B03 Fine structure',fontsize=15,fontweight='bold')
plt.grid()
savefig("B03_fine.png",bbox_inches='tight')
#B03-CO ladder
plt.figure()
plt.loglog(av_b03,l08_b03,label=r'CO(1-0)',linewidth=2)
plt.loglog(av_b03,l09_b03,label=r'CO(2-1)',linewidth=2)
plt.loglog(av_b03,l10_b03,label=r'CO(3-2)',linewidth=2)
plt.loglog(av_b03,l11_b03,label=r'CO(4-3)',linewidth=2)
plt.loglog(av_b03,l12_b03,label=r'CO(5-4)',linewidth=2)
plt.loglog(av_b03,l13_b03,label=r'CO(6-5)',linewidth=2)
plt.loglog(av_b03,l14_b03,label=r'CO(7-6)',linewidth=2)
plt.loglog(av_b03,l15_b03,'--',label=r'CO(8-7)',linewidth=2)
plt.loglog(av_b03,l16_b03,'--',label=r'CO(9-8)',linewidth=2)
plt.loglog(av_b03,l17_b03,'--',label=r'CO(10-9)',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B03 CO ladder',fontsize=15,fontweight='bold')
plt.grid()
savefig("B03_ladder.png",bbox_inches='tight')

#B04-fine structure
plt.figure()
plt.loglog(av_b04,l01_b04,label=r'CII 158$\mu$m',linewidth=2)
plt.loglog(av_b04,l02_b04,label=r'CI 609$\mu$m',linewidth=2)
plt.loglog(av_b04,l04_b04,label=r'CI 320$\mu$m',linewidth=2)
plt.loglog(av_b04,l05_b04,label=r'OI 64$\mu$m',linewidth=2)
plt.loglog(av_b04,l07_b04,label=r'OI 146$\mu$m',linewidth=2)
plt.xlim([1e-6,20])
#plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B04 Fine structure',fontsize=15,fontweight='bold')
plt.grid()
savefig("B04_fine.png",bbox_inches='tight')
#B04-CO ladder
plt.figure()
plt.loglog(av_b04,l08_b04,label=r'CO(1-0)',linewidth=2)
plt.loglog(av_b04,l09_b04,label=r'CO(2-1)',linewidth=2)
plt.loglog(av_b04,l10_b04,label=r'CO(3-2)',linewidth=2)
plt.loglog(av_b04,l11_b04,label=r'CO(4-3)',linewidth=2)
plt.loglog(av_b04,l12_b04,label=r'CO(5-4)',linewidth=2)
plt.loglog(av_b04,l13_b04,label=r'CO(6-5)',linewidth=2)
plt.loglog(av_b04,l14_b04,label=r'CO(7-6)',linewidth=2)
plt.loglog(av_b04,l15_b04,'--',label=r'CO(8-7)',linewidth=2)
plt.loglog(av_b04,l16_b04,'--',label=r'CO(9-8)',linewidth=2)
plt.loglog(av_b04,l17_b04,'--',label=r'CO(10-9)',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel(r'Line emissivity [erg/cm$^3$/s]',fontsize=15)
plt.title('B04 CO ladder',fontsize=15,fontweight='bold')
plt.grid()
savefig("B04_ladder.png",bbox_inches='tight')



#plt.show()
