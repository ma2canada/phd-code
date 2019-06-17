import numpy as np
np.seterr(divide='ignore', invalid='ignore')

def bb(n,K):
    return float(n);
def dd(n,K):
    return float(n*n)/float(K);
def gg(n,K):
    return bb(n,K)+dd(n,K);

def matrixpopulator(i,j,K):
    output=0.0;
    if(i==j):
        #output="-g%i"%i;
        output=-gg(i,K);
    elif(i+1==j):
        #output="b%i"%i;
        output=bb(i,K);
    elif(i==j+1):
        #output="d%i"%i;
        output=dd(i,K);
    return output;

def makethematrix(K):
    #maxn=int(K/2)
    maxn=int(K)
    ourtrix=np.array([[matrixpopulator(i+1,j+1,K) for j in range(maxn)] for i in range(maxn)])
    return ourtrix;

def exitvec1(K):#extinction
    #maxn=int(K/2)
    maxn=int(K)
    outvec=np.array([0. for i in range(maxn)])
    outvec[0]=-dd(1,K)
    #np.insert(outvec,0,-dd(1,K))
    return outvec;
def exitvec2(K):#invasion
    #maxn=int(K/2)
    maxn=int(K)
    outvec=np.array([0. for i in range(maxn)])
    outvec[-1]=-bb(maxn-1,K)
    return outvec;

def probvec1(K):
    return np.linalg.inv(makethematrix(K)).dot(exitvec1(K))
def probvec2(K):#questionable
    #return np.array([1. for i in range(K)])-probvec1(K)
    return np.linalg.inv(makethematrix(K)).dot(exitvec2(K))

def phivec1(K):
    return -np.linalg.inv(makethematrix(K)).dot(probvec1(K))
def phivec2(K):
    return -np.linalg.inv(makethematrix(K)).dot(probvec2(K))

def tauinv(K):
    return phivec2(K)/probvec2(K)
def taufai(K):
    return phivec1(K)/probvec1(K)

#print(dd(1,4))
#print(exitvec1(4))
#temp=np.linalg.inv(makethematrix(4))
#print(temp.dot(makethematrix(4)))
#print (makethematrix(4))
import matplotlib.pyplot as plt
#plt.plot([probvec1(j)[0] for j in range(1,100)],'m--',label='trytry')#tdatumK[-1],
#plt.xlim(0,100);plt.ylim(0,1)
#plt.plot([taufai(j)[0] for j in range(1,100)],'m--',label='trytry')#tdatumK[-1],
plt.plot([tauinv(j)[0] for j in range(1,100)],'m--',label='trytry')#tdatumK[-1],
plt.show()

