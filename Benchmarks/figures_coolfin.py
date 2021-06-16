from pylab import *
import numpy as np
import matplotlib.pyplot as plt


data=np.genfromtxt('B01.cool.fin')
av_b01=data[:,2]
cii_b01=data[:,3]
ci_b01=data[:,4]
oi_b01=data[:,5]
co_b01=data[:,6]

data=np.genfromtxt('B02.cool.fin')
av_b02=data[:,2]
cii_b02=data[:,3]
ci_b02=data[:,4]
oi_b02=data[:,5]
co_b02=data[:,6]

data=np.genfromtxt('B03.cool.fin')
av_b03=data[:,2]
cii_b03=data[:,3]
ci_b03=data[:,4]
oi_b03=data[:,5]
co_b03=data[:,6]

data=np.genfromtxt('B04.cool.fin')
av_b04=data[:,2]
cii_b04=data[:,3]
ci_b04=data[:,4]
oi_b04=data[:,5]
co_b04=data[:,6]

#B01 
plt.figure()
plt.loglog(av_b01,cii_b01,label='CII',linewidth=2)
plt.loglog(av_b01,ci_b01,label='CI',linewidth=2)
plt.loglog(av_b01,oi_b01,label='OI',linewidth=2)
plt.loglog(av_b01,co_b01,label='CO',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Cooling function',fontsize=15)
plt.title('B01',fontsize=15,fontweight='bold')
plt.grid()
savefig("B01_cool.png",bbox_inches='tight')

#B02 
plt.figure()
plt.loglog(av_b02,cii_b02,label='CII',linewidth=2)
plt.loglog(av_b02,ci_b02,label='CI',linewidth=2)
plt.loglog(av_b02,oi_b02,label='OI',linewidth=2)
plt.loglog(av_b02,co_b02,label='CO',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Cooling function',fontsize=15)
plt.title('B02',fontsize=15,fontweight='bold')
plt.grid()
savefig("B02_cool.png",bbox_inches='tight')

#B03 
plt.figure()
plt.loglog(av_b03,cii_b03,label='CII',linewidth=2)
plt.loglog(av_b03,ci_b03,label='CI',linewidth=2)
plt.loglog(av_b03,oi_b03,label='OI',linewidth=2)
plt.loglog(av_b03,co_b03,label='CO',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Cooling function',fontsize=15)
plt.title('B03',fontsize=15,fontweight='bold')
plt.grid()
savefig("B03_cool.png",bbox_inches='tight')

#B04 
plt.figure()
plt.loglog(av_b04,cii_b04,label='CII',linewidth=2)
plt.loglog(av_b04,ci_b04,label='CI',linewidth=2)
plt.loglog(av_b04,oi_b04,label='OI',linewidth=2)
plt.loglog(av_b04,co_b04,label='CO',linewidth=2)
plt.xlim([1e-6,20])
plt.legend(loc='best')
plt.xlabel('Av [mag]',fontsize=15)
plt.ylabel('Cooling function',fontsize=15)
plt.title('B04',fontsize=15,fontweight='bold')
plt.grid()
savefig("B04_cool.png",bbox_inches='tight')


#plt.show()


