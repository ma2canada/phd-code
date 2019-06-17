import numpy as np
import matplotlib.pyplot as plt

chooseK=50.; choosedelta=9.5; chooseq=.01;

def gg(n,q,d,KK):
    return (1.+d)*n - 1.*q*n*n/KK + ( d*n + (1.-q)*n*n/KK )
#    return (1.+2.*d)*n + (1.-2.*q)*n*n/KK

def asdarg(n,q,d,KK):
    return (1. - 2.*q)*n + KK*(1. + 2.*d);
def potential(n,q,d,KK):
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

nrange = np.arange(.2,1.1*chooseK,.2);
''' these show that the arguments of the log do go negative, but only some distance after K
dumbfig1=plt.figure(1)
plt.plot(nrange,asdarg(nrange,chooseq,choosedelta,chooseK))
plt.show()
dumbfig2=plt.figure(2)
plt.plot(nrange,gg(nrange,chooseq,choosedelta,chooseK))
plt.show()
'''
figgypud=plt.figure(3)
ax = figgypud.add_subplot(111)#,aspect='equal')
plt.plot(nrange,potential(nrange,chooseq,choosedelta,chooseK))
#plt.plot(findMaxAndMin(chooseq,choosedelta,chooseK),[potential(findMaxAndMin(chooseq,choosedelta,chooseK)[0],chooseq,choosedelta,chooseK),potential(findMaxAndMin(chooseq,choosedelta,chooseK)[1],chooseq,choosedelta,chooseK)],'o',color='red')
plt.plot([findMaxAndMin(chooseq,choosedelta,chooseK)[0]],[potential(findMaxAndMin(chooseq,choosedelta,chooseK)[0],chooseq,choosedelta,chooseK)],'o',color='red')
plt.plot([findMaxAndMin(chooseq,choosedelta,chooseK)[1]],[potential(findMaxAndMin(chooseq,choosedelta,chooseK)[1],chooseq,choosedelta,chooseK)],'o',color='green')
''' red dot should be higher than green dot '''
plt.ylabel('potential',fontsize=14)
plt.xlabel('$n$',fontsize=14)
#forceAspect(ax,aspect=1)
#ax.set_aspect('equal')
plt.show()
print(depth(chooseq,choosedelta,chooseK))

#print(findMaxAndMin(chooseq,choosedelta,chooseK))
#print([potential(findMaxAndMin(chooseq,choosedelta,chooseK)[0],chooseq,choosedelta,chooseK),potential(findMaxAndMin(chooseq,choosedelta,chooseK)[1],chooseq,choosedelta,chooseK)])

#cc=np.array([1,2,3]); bb=np.array([5])
#print(np.concatenate((cc,bb)))
#a = np.array([[1, 2], [3, 4]])
#b = np.array([[5, 6]])
#print(np.concatenate((a, b), axis=0))
