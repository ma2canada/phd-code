import matplotlib.pyplot as plt
import numpy as np
#import pylab as pl
from matplotlib import rc
rc('text', usetex=True)

gdata = np.arange(0.01,1.-0.01,0.001)
gdatalo = np.arange(0.01,0.5+0.002,0.001)
gdatahi = np.arange(0.5-0.002,1.0-0.01,0.001)
xaxislo = np.array([1E-1 for i in range(len(gdatalo))])
xaxishi = np.array([1E-1 for i in range(len(gdatahi))])
toplo = np.array([1E2 for i in range(len(gdatalo))])
tophi = np.array([1E2 for i in range(len(gdatahi))])
y2datlo = 1./(1.-gdatalo)
y3datlo = 1./gdatalo
y2dathi = 1./(1.-gdatahi)
y3dathi = 1./gdatahi

c1='yellow'#(0.3,0.4,0.1)#'#842F32'
c2='chartreuse'#(0.0,0.4,0.5)#'#878092'
c3='cornflowerblue'#(0.3,0.0,0.4)#'#E3E3CD'
#plt.figure()
#plt.Axes.fill_between(xdata,y1dat,y2dat)
regionfig,(ax) = plt.subplots()
ax.fill_between(gdatalo,xaxislo,y2datlo,facecolor=c1)#Nv<1/g
ax.fill_between(gdatalo,toplo,y3datlo,facecolor=c2)#Nv>1/g
ax.fill_between(gdatahi,xaxishi,y3dathi,facecolor=c1)#Nv<1/g
ax.fill_between(gdatahi,tophi,y2dathi,facecolor=c2)#Nv>1/1-g
ax.fill_between(gdatalo,y2datlo,y3datlo,facecolor=c3)#intermediate
ax.fill_between(gdatahi,y2dathi,y3dathi,facecolor=c3)#intermediate
fntsz=16
#ax.axvline(x=0, color='k',lw=0.5)
#ax.axhline(y=0, color='k',lw=0.5)
ax.text(0.15, 10.**1.2, 'population maintenance', fontsize=fntsz)#, ha='center'
ax.text(0.05, 10.**0.3, 'intermediate', fontsize=fntsz)
ax.text(0.15, 10.**-0.3, 'fixated or extinct', fontsize=fntsz)
ax.text(0.25, 10.**1.7, r'$\uparrow$ frequent immigration', fontsize=fntsz)
ax.text(0.30, 10.**-0.8, r'$\downarrow$ rare immigration', fontsize=fntsz)
plt.ylabel(r'$N\nu$',fontsize=fntsz+4)
plt.xlabel(r'$g$',fontsize=fntsz+4)
plt.yscale("log")
plt.ylim(1E-1,1E2)
plt.xlim(0,1)
#ax.set_aspect(aspect='equal')
#plt.legend()
plt.show()
regionfig.savefig("ch3regimes.pdf",bbox_inches='tight')
