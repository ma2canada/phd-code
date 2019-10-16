#1D logistic seeing how eigenvalues vary with K
# eigvals failing me at tau~1e11, regardless of the shift I put on it - not a stack overflow?!?

from scipy.sparse.linalg import eigs
import datetime
import scipy.sparse as sparse
from numpy.linalg import eigvals
import numpy as np

kconst = 1;
buffr = 5;
maxi = kconst*buffr;

def bn(n,k):
    return np.float_(n);
def dn(n,k):
    return np.float_(n*float(n)/k);
def gg(n,k):
    return np.float_(bn(n,k)+dn(n,k));

def makesparsetrix1D(kc,maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi):
      for j in range(maxi):
            if(i == j and i == maxi-1):
                rowvec.append(i);
                colvec.append(j);
                datvec.append(-dn(j + 1,kc));#=modified g for reflecting
            elif(i == j):
                rowvec.append(i);
                colvec.append(j);
                datvec.append(-gg(j + 1,kc));
            elif(i-1 == j):
                rowvec.append(i);
                colvec.append(j);
                datvec.append(bn(j + 1,kc));
            elif(i+1 == j):
                rowvec.append(i);
                colvec.append(j);
                datvec.append(dn(j + 1,kc));
    trix = sparse.coo_matrix((datvec,(rowvec,colvec)), shape=(maxi,maxi))
    return trix

def makedensetrix1D(kc,maxi):
    mat = [[0. for j in range(maxi)] for i in range(maxi)];
    for i in range(maxi):
      for j in range(maxi):
            if(i == j and i == maxi-1):
                mat[i][j]=-dn(j + 1,kc)*(1.e0);
            elif(i == j):
                mat[i][j]=-gg(j + 1,kc)*(1.e0);
            elif(i-1 == j):
                mat[i][j]=bn(j + 1,kc)*(1.e0);
            elif(i+1 == j):
                mat[i][j]=dn(j + 1,kc)*(1.e0);
    return mat

def guessvec(kc,maxi):
    temp = np.zeros(maxi); temp[kc-1]=0.3;
    return temp

fille = open("sparseevalvK0.txt",'w'); #fille.write("K\tsmallest eval\tsecond smallest eval\n");
#sometimes just... stops? without error messages, just doesn't complete the loop
for kount in range(23):
    kconst=kount*4+4;
    temp3=makesparsetrix1D(kconst,kconst*buffr)#*(1.e+10)
##    temp3=makedensetrix1D(kconst,kconst*buffr)
    print "matrix %i acquired\t"%kconst+str(datetime.datetime.now())
    smallev = eigs(temp3,k=2,which='SM',return_eigenvectors=False,maxiter=10000,tol=1E-5);#,tol=1E-2
##    smallev, smallevec = eigs(temp3,k=2,which='SM',return_eigenvectors=True);#,v0=guessvec(kconst,kconst*buffr)
##    smallev, smallevec = eigs(temp3,k=2,sigma=-0.00001,which='LM',return_eigenvectors=True);#NOT YET SUPPORTED!!
##    smallev = eigvals(temp3.todense()); smallev=sorted(smallev);
##    smallev = eigvals(temp3); smallev=sorted(smallev);
    print smallev.real
##    print smallevec[:,1].real
##    print max(smallevec[:,1].real)
    boop=open("evaldistr/K%i.txt"%kconst,'w');
    for itm in smallev:
        boop.write(str(itm.real)+"\t");
    boop.close();
    fille.write("%i\t%e\t%e\n" %(kconst,-1/smallev[-1].real,-1/smallev[-2].real));#order reverses when eigenvecs calc'ed
fille.close();
print datetime.datetime.now();
