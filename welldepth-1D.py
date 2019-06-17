import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpmath import hyp3f2

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

chooseK=50.;

def gg(n,q,d,KK):
    return (1.+d)*n - 1.*q*n*n/KK + ( d*n + (1.-q)*n*n/KK )
#    return (1.+2.*d)*n + (1.-2.*q)*n*n/KK

def potential(n,q,d,KK):
    '''don't need to worry if it complains about bad log values, this is only for n>K and is irrelevant for Kramers '''
#def potential(n,q,londie,KK):
#    d=10.**londie;
#    if((1.+2.*d)*KK + (1.-2.*q)*n>0.):
#    if(True):
    asdf=-4.*KK*(1. + d - q)*np.log((1. - 2.*q)*n + KK*(1. + 2.*d))/(1. - 2.*q)**2 + 2.*n/(1. - 2.*q) + np.log(gg(n,q,d,KK));
#    elif((1.+2.*d)*KK + (1.-2.*q)*n<0.):
#        asdf=-4.*KK*(1. + d - q)*np.log(-(1. - 2.*q)*n - KK*(1. + 2.*d))/(1. - 2.*q)**2 + 2.*n/(1. - 2.*q) + np.log(gg(n,q,d,KK))
#    else:
#        print("bad bad not good");
#        asdf=0;
    return asdf
def findMaxAndMin(q,d,KK):
    resolution=.2
    nrange = np.arange(resolution,2*KK,resolution)
    potrange = potential(nrange,q,d,KK)
    #find the first differences
    potrange=np.concatenate((potrange,[potrange[-1]]))-np.concatenate(([potrange[0]-1.],potrange))
    maxposition=-1;minposition=0;
    for i,elem in enumerate(potrange):
        if(elem<0. and maxposition==-1):
            maxposition=i-1;
        if(maxposition>=0 and elem>0.):
            minposition=i-1;
            break;#okay but what happens if it's a monotonic function?
    if(maxposition==0):
        maxposition=1;
    if(minposition==0):
        minposition=int(KK);
    return [nrange[maxposition],nrange[minposition]]
def depth(q,d,KK):
    [top,btm]=findMaxAndMin(q,d,KK)
    #bottom=potential(KK,q,d,KK)
    #top=potential(1.,q,d,KK)
    #return top-bottom
    return potential(top,q,d,KK)-potential(btm,q,d,KK)

num = 20; dx = 1./num;
#ddd = np.logspace(0.1, 10., num)
lnd = np.arange(-1.0, 1.0, 2*dx)#the problem comes for q=1/2 exactly
cue = np.arange(-0.0001, 1.0-0.0001, dx)
X, Y = np.meshgrid(cue,lnd)
stochasticity = np.linspace(0.01+.0001, .99+.0001, 100)
variability = np.logspace(-1.0, 1.0, 1000)
X,Y=np.meshgrid(stochasticity,variability)
#Z = depth(X, 10.**Y, chooseK)

"""
''' not great "correct" solution with which to compare '''
def mte1Dsum_tau1(stoch,delta,cap):
    return abs(float(2*cap*hyp3f2(1,1,1-cap/stoch-delta*cap/(2*stoch),2,(2+delta*cap/2-2*stoch)/(stoch-2),stoch/(stoch-1))/(2+delta*cap-2*stoch)))
Z2 = Z;
#for i1,x1 in enumerate(X):
#    for i2,x2 in enumerate(Y):
for i1 in range(len(X)):
    for i2 in range(len(X[0])):
        Z2[i1][i2]=mte1Dsum_tau1(X[i1][i2], 10.**Y[i1][i2], chooseK)
"""
''' hopefully better solution with which to compare '''
MTE2 = np.load("C:/Users/lenov/Documents/GitHub/shift-stochasticity/data/heat_MTE_K100_log.npy")
#stochasticity = np.linspace(0.01, .99, 100)
#variability = np.logspace(-1.0, 1.0, 1000)
Z2=np.asarray(MTE2)
Z = [[0. for j1 in Z2[0]] for j2 in Z2];
'''print(len(np.asarray(MTE2)))
print(len(np.asarray(MTE2)[0]))
print(len(stochasticity))
print(len(variability))
print(len(Z))
print(len(Z[0]))'''
for i1,elem1 in enumerate(stochasticity):
    for i2,elem2 in enumerate(variability):
        Z[i2][i1]=depth(elem2,10.**elem1,chooseK)

#fig, ax  = plt.subplots()
plt.figure()
plt.xlabel('$q$',fontsize=16)
plt.ylabel('$\ln(\delta)$',fontsize=16)
minum=-.0000; maxum=.0014; numum=7.; stepum = (maxum-minum)/numum;
levels = list(np.arange(minum,maxum,stepum))#as informed by looking at the graph without setting levels
CS = plt.contourf(X, Y, Z/Z2,levels,cmap=plt.cm.YlGnBu)#
#plt.clabel(CS, inline=1, fontsize=10)
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel("ratio of Kramer's depth to true extinction time")
#CS.ax.set_yscale("log")
plt.show()

minval=1.E16;
for elem in Z2:
    minval=min(minval,min(elem))
print(minval,"=min(Z2) is the minimum of Z2, by which we divide Z,")
#print("and this is Z2 itself: ", Z2)
minvalZ=1.E16;
for elem in Z:
    minvalZ=min(minvalZ,min(elem))
print(minvalZ,"=min(Z) is the minimum of Z, which we divide by Z2,")
#print("and this is Z2 itself: ", Z2)

print(Z[-1][0:8])
print(Z2[-1][0:8])
print((Z/Z2)[-1][0:8])
