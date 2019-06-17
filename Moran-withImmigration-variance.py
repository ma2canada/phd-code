import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
import numpy as np

def varinf(N,nu,g):
    return g*(1-g)*N*N/(1+nu*(N-1))

def vart(N,nu,g,t,mu0,var0):
    A=(1+g*nu-g*(1-nu)/N)*N*N*(mu0-g*N)/(N*nu+2-2*nu)
    B=(g*N-mu0)**2
    C=var0-varinf(N,nu,g)+B+(g*N-mu0)*(2-nu)*(1-2*g)/(N*nu+2-2*nu)
    return varinf(N,nu,g)+A*np.exp(-nu*t/N)-B*np.exp(-2*nu*t/N)+C*np.exp(-2*(nu+(1-nu)/N)*t/N)

currentN=100;
currentg=0.1;
currentNnu=8.;
currentnu=float(currentNnu)/float(currentN);
currentmu0=currentg*currentN*1.1;
currentvar0=0.01*currentN**2;

'''
t=np.linspace(0,12000,1000)

plt.figure(1)
#plt.semilogy(t,vart(currentN,currentnu,currentg,t,currentmu0,currentvar0))#/currentN**2
plt.semilogy(t,vart(currentN,currentnu,currentg,t,currentmu0,currentvar0))#/currentN**2
plt.semilogy(t,np.exp(-currentnu*t/currentN)*(vart(currentN,currentnu,currentg,0.,currentmu0,currentvar0)-varinf(currentN,currentnu,currentg))+varinf(currentN,currentnu,currentg))
plt.xlabel('time $t$')
plt.ylabel('variance')
#plt.legend()
#plt.savefig('MoranVarianceTimed.pdf', bbox_inches='tight')
plt.show()
'''

gspace=np.linspace(0,1,200)
#nuspace=np.logspace(-2,-0,num=100)
#nuspace=np.geomspace(0.01,1.0,10+1)
nuspace=np.geomspace(0.0005,0.1,200+1)
#IHateArrays=[[varinf(currentN,nu,g) for g in gspace] for nu in nuspace]
IHateArrays=[[varinf(currentN,nu,g)/currentN**2 for g in gspace] for nu in nuspace]

#plt.figure(2)
im=plt.contourf(gspace,nuspace,IHateArrays,cmap='Greens',levels=100)
plt.xlabel('metapopulation abundance, $g$')
plt.ylabel(r'immigration rate, $\nu$')
plt.yscale('log')
'''
firsthalf=np.linspace(0.,0.5,100)
firsthalfshort=np.linspace(0.1,0.5,100)
seconhalf=np.linspace(0.5,1.,100)
seconhalfshort=np.linspace(0.5,0.9,100)
firsthalfY=1./firsthalfshort/currentN
seconhalfY=1./seconhalf/currentN
firsthalfW=1./(1.-firsthalf)/currentN
seconhalfW=1./(1.-seconhalfshort)/currentN
plt.plot(firsthalfshort,firsthalfY)
plt.plot(seconhalf,seconhalfY)
plt.plot(firsthalf,firsthalfW)
plt.plot(seconhalfshort,seconhalfW)
'''
ax = plt.gca() #gca = get current axis
# create colorbar
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel(r'normalized variance',va="top",fontsize=12)#,rotation=-90,labelpad=20
cbar.ax.tick_params(labelsize=12)

#plt.savefig("MoranVariance.pdf",bbox_inches='tight')
plt.show()

