#generates the transition matrix associated with a process (SIR w/ vital)
import birthdeathrates2020 as bdr
##from scipy.linalg import expm
from numpy import dot,ceil
import scipy.sparse as sparse

#n is Sus, m is Inf
'''one indexing'''
def nfromindex(ind,maxi):
    return ind%maxi+1
def mfromindex(ind,maxi):
    return ind//maxi+1
def indexfromnm(n,m,maxi):
    return (m-1)*maxi+(n-1)
'''zero indexing'''
def nfromindex0(ind,maxi):
    return ind%(maxi+1)
def mfromindex0(ind,maxi):
    return ind//(maxi+1)
def index0fromnm(n,m,maxi):
    return m*(maxi+1)+n
'''S from zero, I from one'''
def nfromindexJ(ind,maxi):
    return ind%(maxi+1)
def mfromindexJ(ind,maxi):
    return ind//(maxi+1)+1
def indexJfromnm(n,m,maxi):
    #print((m-1)*(maxi+1)+n)
    return (m-1)*(maxi+1)+n

def maketrix(a,kc,maxi):
    trix = [[0 for y in range(maxi*maxi)] for x in range(maxi*maxi)]
    for i in range(maxi):
      for j in range(maxi):
        for k in range(maxi):
          for l in range(maxi):
            if(i == j and k == l and j + 1 == maxi and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gnm(l+1, j + 1,a,kc);
            elif(i == j and k == l and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gn(l+1, j + 1,a,kc);
            elif(i == j and k == l and j+1 == maxi):
                trix[maxi*i + k][maxi*j + l]=-bdr.gm(l+1, j + 1,a,kc);
            elif(i == j and k == l):
                trix[maxi*i + k][maxi*j + l]=-bdr.gg(l+1, j + 1,a,kc);
            elif(i == j and k-1 == l):
                trix[maxi*i + k][maxi*j + l]= bdr.bn(l+1, j + 1,a,kc);
            elif(i == j and k+1 == l):
                trix[maxi*i + k][maxi*j + l]= bdr.dn(l+1, j + 1,a,kc);
            elif(i-1 == j and k == l):
                trix[maxi*i + k][maxi*j + l]= bdr.bm(l+1, j + 1,a,kc);
            elif(i+1 == j and k == l):
                trix[maxi*i + k][maxi*j + l]= bdr.dm(l+1, j + 1,a,kc);
    return trix

def pinit(a,kc,maxi):
    n0=int(round(kc/(1.+a)));
    initcond=((n0 - 1)*maxi + n0-1)
    pinit = [0 for y in range(maxi*maxi)]; pinit[initcond]=1.;
    return pinit

def maketrix_test(params,maxi):
    trix = [[' 0' for y in range(maxi*maxi)] for x in range(maxi*maxi)]
    for i in range(maxi):
      for j in range(maxi):
        for k in range(maxi):
          for l in range(maxi):
            if(i == j and k == l and j + 1 == maxi and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]='nm';
            elif(i == j and k == l and l+1 == maxi):
                trix[maxi*i + k][maxi*j + l]='gn';
            elif(i == j and k == l and j+1 == maxi):
                trix[maxi*i + k][maxi*j + l]='gm';
            elif(i == j and k == l):
                trix[maxi*i + k][maxi*j + l]='gg';
            elif(i == j and k-1 == l):
                trix[maxi*i + k][maxi*j + l]='bn';
            elif(i == j and k+1 == l):
                trix[maxi*i + k][maxi*j + l]='dn';
            elif(i-1 == j and k+1 == l):#here's the change
                trix[maxi*i + k][maxi*j + l]='bm';
            elif(i-1 == j and k == l):#here's unchanged
                trix[maxi*i + k][maxi*j + l]='BM';
            elif(i+1 == j and k == l):
                trix[maxi*i + k][maxi*j + l]='dm';
    return trix

##def pfinal(a,kc,maxi,time):
##    matr=maketrix(a,kc,maxi);
##    pfinal = [0 for y in range(maxi*maxi)];
##    pfinal = dot( expm([[matr[x][y]*time for x in range(maxi*maxi)] for y in range(maxi*maxi)]).T, pinit(a,kc,maxi) );
##    # I think when multiplying by the time this creates a whole new matrix
##    return pfinal

#this is most definitely NOT the best way to make it, it could be one loop or five
#NTS: check out sparse.diags
def makesparsetrix(params,maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi):
      for j in range(maxi):
        for k in range(maxi):
          for l in range(maxi):
            if(i == j and k == l and j + 1 == maxi and l+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gnm(l+1, j + 1,params));
            elif(i == j and k == l and l+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gn(l+1, j + 1,params));
            elif(i == j and k == l and j+1 == maxi):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gm(l+1, j + 1,params));
            elif(i == j and k == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(-bdr.gg(l+1, j + 1,params));
            elif(i == j and k-1 == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.bn(l+1, j + 1,params));
            elif(i == j and k+1 == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.dn(l+1, j + 1,params));
            elif(i-1 == j and k+1 == l):                   #I THINK this (k->k+1) is the right change to account for birth AND death
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.bm(l+1, j + 1,params));
            elif(i+1 == j and k == l):
                rowvec.append(maxi*i + k);
                colvec.append(maxi*j + l);
                datvec.append(+bdr.dm(l+1, j + 1,params));
    trix = sparse.coo_matrix((datvec,(rowvec,colvec)), shape=(maxi*maxi,maxi*maxi))
    return trix
'''
#using sparse.diags
def makesmartrix(params,maxi):
    diagmain=[];
    for i in range(maxi*maxi):
        diagmain.append(-bdr.gg(nfromindex(i,maxi),mfromindex(i,maxi),params))
    for j in range(maxi):
        diagmain[indexfromnm(j+1,maxi,maxi)]=-bdr.gm(j+1,maxi,params)
        diagmain[indexfromnm(maxi,j+1,maxi)]=-bdr.gn(maxi,j+1,params)
    diagmain[indexfromnm(maxi,maxi,maxi)]=-bdr.gnm(maxi,maxi,params)
    diagbn=[0 for i in range(maxi*maxi-1)]
        
    
        trix = sparse.diags([diagmain,diagbn,diagbm,diagdn,diagdm],[0,+1,+maxi,-1,-maxi])
    return trix
'''

def makeJemtrix(params,maxi):
    rowvec = []; colvec = []; datvec = [];
    '''print("(0,1)")'''
    #none one = (S,I) = (0,1)
    rowvec.append(indexJfromnm(0,1,maxi));
    colvec.append(indexJfromnm(0,1,maxi));
    datvec.append(-bdr.gg(0,1,params));
    rowvec.append(indexJfromnm(0,1,maxi));
    colvec.append(indexJfromnm(0+1,1,maxi));
    datvec.append(bdr.dn(1,1,params));
    rowvec.append(indexJfromnm(0,1,maxi));
    colvec.append(indexJfromnm(0,1+1,maxi));
    datvec.append(bdr.dm(0,1+1,params));
    #rowvec.append(indexJfromnm(1,0,maxi));
    #colvec.append(indexJfromnm(1,0,maxi));
    #datvec.append(bdr.bm(1,0,params));#state out of bounds
    '''print("(1,1)")'''
    #one infection, one susc
    rowvec.append(indexJfromnm(1,1,maxi));
    colvec.append(indexJfromnm(1,1,maxi));
    datvec.append(-bdr.gg(1,1,params));
    rowvec.append(indexJfromnm(1,1,maxi));
    colvec.append(indexJfromnm(1+1,1,maxi));
    datvec.append(bdr.dn(1+1,1,params));
    rowvec.append(indexJfromnm(1,1,maxi));
    colvec.append(indexJfromnm(1,1+1,maxi));
    datvec.append(bdr.dm(1,1+1,params));
    rowvec.append(indexJfromnm(1,1,maxi));
    colvec.append(indexJfromnm(1-1,1,maxi));
    datvec.append(bdr.bn(1-1,1,params));
    for j in range(2,(maxi+1)-1):
        '''print("(S,1)")'''
        #one infection
        rowvec.append(indexJfromnm(j,1,maxi));
        colvec.append(indexJfromnm(j,1,maxi));
        datvec.append(-bdr.gg(j,1,params));
        rowvec.append(indexJfromnm(j,1,maxi));
        colvec.append(indexJfromnm(j+1,1,maxi));
        datvec.append(bdr.dn(j+1,1,params));
        rowvec.append(indexJfromnm(j,1,maxi));
        colvec.append(indexJfromnm(j,1+1,maxi));
        datvec.append(bdr.dm(j,1+1,params));
        #rowvec.append(index0fromnm(i,1,maxi));
        #colvec.append(index0fromnm(i,1,maxi)-1);
        #datvec.append(bdr.bn(i+1,1-1,params));
        rowvec.append(indexJfromnm(j,1,maxi));
        colvec.append(indexJfromnm(j-1,1,maxi));
        datvec.append(bdr.bn(j-1,1,params));
        '''print("(0,I)")'''
        #no susceptible!!!!
        rowvec.append(indexJfromnm(0,j,maxi));
        colvec.append(indexJfromnm(0,j,maxi));
        datvec.append(-bdr.gg(0,j,params));
        rowvec.append(indexJfromnm(0,j,maxi));
        colvec.append(indexJfromnm(0+1,j,maxi));
        datvec.append(bdr.dn(1,j,params));
        rowvec.append(indexJfromnm(0,j,maxi));
        colvec.append(indexJfromnm(0,j+1,maxi));
        datvec.append(bdr.dm(0,j+1,params));
        rowvec.append(indexJfromnm(0,j,maxi));
        colvec.append(indexJfromnm(0+1,j-1,maxi))#-(maxi+1)+1);
        datvec.append(bdr.bm(0+1,j-1,params));
        for i in range(1,maxi+1-1):
            '''print("(S,I)")'''
            rowvec.append(indexJfromnm(i,j,maxi));
            colvec.append(indexJfromnm(i,j,maxi));
            datvec.append(-bdr.gg(i,j,params));
            rowvec.append(indexJfromnm(i,j,maxi));
            colvec.append(indexJfromnm(i+1,j,maxi));
            datvec.append(bdr.dn(i+1,j,params));
            rowvec.append(indexJfromnm(i,j,maxi));
            colvec.append(indexJfromnm(i,j+1,maxi))#+(maxi+1));
            datvec.append(bdr.dm(i,j+1,params));
            rowvec.append(indexJfromnm(i,j,maxi));
            colvec.append(indexJfromnm(i-1,j,maxi));
            datvec.append(bdr.bn(i-1,j,params));
            rowvec.append(indexJfromnm(i,j,maxi));
            colvec.append(indexJfromnm(i+1,j-1,maxi))#-(maxi+1)+1);
            datvec.append(bdr.bm(i+1,j-1,params));
        '''print("(S,N)")'''
        #full infection
        rowvec.append(indexJfromnm(j,maxi,maxi));
        colvec.append(indexJfromnm(j,maxi,maxi));
        datvec.append(-bdr.gm(j,maxi,params));
        rowvec.append(indexJfromnm(j,maxi,maxi));
        colvec.append(indexJfromnm(j+1,maxi,maxi));
        datvec.append(bdr.dn(j+1,maxi,params));
        rowvec.append(indexJfromnm(j,maxi,maxi));
        colvec.append(indexJfromnm(j-1,maxi,maxi));
        datvec.append(bdr.bn(j-1,maxi,params));
        rowvec.append(indexJfromnm(j,maxi,maxi));
        colvec.append(indexJfromnm(j+1,maxi-1,maxi))#-(maxi+1)+1);
        datvec.append(bdr.bm(j+1,maxi-1,params));
        '''print("(N,I)")'''
        #full susceptible
        rowvec.append(indexJfromnm(maxi,j,maxi));
        colvec.append(indexJfromnm(maxi,j,maxi));
        datvec.append(-bdr.gn(maxi,j,params));
        rowvec.append(indexJfromnm(maxi,j,maxi));
        colvec.append(indexJfromnm(maxi,j+1,maxi))#+(maxi+1));
        datvec.append(bdr.dm(maxi,j+1,params));
        rowvec.append(indexJfromnm(maxi,j,maxi));
        colvec.append(indexJfromnm(maxi-1,j,maxi));
        datvec.append(bdr.bn(maxi-1,j,params));
        #rowvec.append(index0fromnm(maxi,i,maxi));
        #colvec.append(index0fromnm(maxi,i,maxi)-(maxi+1)+1);
        #datvec.append(bdr.bm(maxi+1,i-1,params));
    '''print("(0,N)")'''
    #full infection, no susc
    rowvec.append(indexJfromnm(0,maxi,maxi));
    colvec.append(indexJfromnm(0,maxi,maxi));
    datvec.append(-bdr.gm(0,maxi,params));
    rowvec.append(indexJfromnm(0,maxi,maxi));
    colvec.append(indexJfromnm(0+1,maxi,maxi));
    datvec.append(bdr.dn(0+1,maxi,params));
    rowvec.append(indexJfromnm(0,maxi,maxi));
    colvec.append(indexJfromnm(0+1,maxi-1,maxi))#-(maxi+1)+1);
    datvec.append(bdr.bm(0+1,maxi-1,params));
    '''print("(1,N)")'''
    #full infection, one susc
    rowvec.append(indexJfromnm(1,maxi,maxi));
    colvec.append(indexJfromnm(1,maxi,maxi));
    datvec.append(-bdr.gm(1,maxi,params));
    rowvec.append(indexJfromnm(1,maxi,maxi));
    colvec.append(indexJfromnm(1+1,maxi,maxi));
    datvec.append(bdr.dn(1+1,maxi,params));
    rowvec.append(indexJfromnm(1,maxi,maxi));
    colvec.append(indexJfromnm(1-1,maxi,maxi));
    datvec.append(bdr.bn(1-1,maxi,params));
    rowvec.append(indexJfromnm(1,maxi,maxi));
    colvec.append(indexJfromnm(1+1,maxi-1,maxi))#-(maxi+1)+1);
    datvec.append(bdr.bm(1+1,maxi-1,params));
    '''print("(N,1)")'''
    #full susceptible, one infe
    rowvec.append(indexJfromnm(maxi,1,maxi));
    colvec.append(indexJfromnm(maxi,1,maxi));
    datvec.append(-bdr.gn(maxi,j,params));
    rowvec.append(indexJfromnm(maxi,1,maxi));
    colvec.append(indexJfromnm(maxi,1+1,maxi))#+(maxi+1));
    datvec.append(bdr.dm(maxi,1+1,params));
    rowvec.append(indexJfromnm(maxi,1,maxi));
    colvec.append(indexJfromnm(maxi-1,1,maxi));
    datvec.append(bdr.bn(maxi-1,1,params));
    '''print("(N,N)")'''
    #full full
    rowvec.append(indexJfromnm(maxi,maxi,maxi));
    colvec.append(indexJfromnm(maxi,maxi,maxi));
    datvec.append(-bdr.gnm(maxi,maxi,params));
    rowvec.append(indexJfromnm(maxi,maxi,maxi));
    colvec.append(indexJfromnm(maxi-1,maxi,maxi));
    datvec.append(bdr.bn(maxi-1,maxi,params));
    '''print("out of",(maxi)*(maxi+1))'''
    trix = sparse.coo_matrix((datvec,(rowvec,colvec)), shape=((maxi)*(maxi+1),(maxi)*(maxi+1)))
    #print("dummy, the (square) shape is %i"%(maxi*maxi+maxi))
    return trix

def make0emtrix(params,maxi):
    rowvec = []; colvec = []; datvec = [];
    #none none
    rowvec.append(index0fromnm(0,0,maxi));
    colvec.append(index0fromnm(0,0,maxi));
    datvec.append(-bdr.gg(0,0,params));
    rowvec.append(index0fromnm(0,0,maxi));
    colvec.append(index0fromnm(0,0,maxi)+1);#or index0fromnm(0+1,0,maxi)
    datvec.append(bdr.dn(1,0,params));
    rowvec.append(index0fromnm(0,0,maxi));
    colvec.append(index0fromnm(0,0+1,maxi));
    datvec.append(bdr.dm(0,1,params));
    for i in range(1,maxi+1-1):
        #no infection
        rowvec.append(index0fromnm(i,0,maxi));
        colvec.append(index0fromnm(i,0,maxi));
        datvec.append(-bdr.gg(i,0,params));
        rowvec.append(index0fromnm(i,0,maxi));
        colvec.append(index0fromnm(i,0,maxi)+1);
        datvec.append(bdr.dn(i+1,0,params));
        rowvec.append(index0fromnm(i,0,maxi));
        colvec.append(index0fromnm(i,0,maxi)+(maxi+1));
        datvec.append(bdr.dm(i,1,params));
        rowvec.append(index0fromnm(i,0,maxi));
        colvec.append(index0fromnm(i,0,maxi)-1);
        datvec.append(bdr.bn(i-1,0,params));
        #no susceptible
        rowvec.append(index0fromnm(0,i,maxi));
        colvec.append(index0fromnm(0,i,maxi));
        datvec.append(-bdr.gg(0,i,params));
        rowvec.append(index0fromnm(0,i,maxi));
        colvec.append(index0fromnm(0,i,maxi)+1);
        datvec.append(bdr.dn(1,i,params));
        rowvec.append(index0fromnm(0,i,maxi));
        colvec.append(index0fromnm(0,i,maxi)+(maxi+1));
        datvec.append(bdr.dm(0,i+1,params));
        rowvec.append(index0fromnm(0,i,maxi));
        colvec.append(index0fromnm(0+1,i-1,maxi))#-(maxi+1)+1);
        datvec.append(bdr.bm(0+1,i-1,params));
        for j in range(1,(maxi+1)-1):
            rowvec.append(index0fromnm(i,j,maxi));
            colvec.append(index0fromnm(i,j,maxi));
            datvec.append(-bdr.gg(i,j,params));
            rowvec.append(index0fromnm(i,j,maxi));
            colvec.append(index0fromnm(i,j,maxi)+1);
            datvec.append(bdr.dn(i+1,j,params));
            rowvec.append(index0fromnm(i,j,maxi));
            colvec.append(index0fromnm(i,j+1,maxi))#+(maxi+1));
            datvec.append(bdr.dm(i,j+1,params));
            rowvec.append(index0fromnm(i,j,maxi));
            colvec.append(index0fromnm(i,j,maxi)-1);
            datvec.append(bdr.bn(i-1,j,params));
            rowvec.append(index0fromnm(i,j,maxi));
            colvec.append(index0fromnm(i+1,j-1,maxi))#-(maxi+1)+1);
            datvec.append(bdr.bm(i+1,j-1,params));
        #full infection
        rowvec.append(index0fromnm(i,maxi,maxi));
        colvec.append(index0fromnm(i,maxi,maxi));
        datvec.append(-bdr.gm(i,maxi,params));
        rowvec.append(index0fromnm(i,maxi,maxi));
        colvec.append(index0fromnm(i,maxi,maxi)+1);
        datvec.append(bdr.dn(i+1,maxi,params));
        rowvec.append(index0fromnm(i,maxi,maxi));
        colvec.append(index0fromnm(i,maxi,maxi)-1);
        datvec.append(bdr.bn(i-1,maxi,params));
        rowvec.append(index0fromnm(i,maxi,maxi));
        colvec.append(index0fromnm(i+1,maxi-1,maxi))#-(maxi+1)+1);
        datvec.append(bdr.bm(i+1,maxi-1,params));
        #full susceptible
        rowvec.append(index0fromnm(maxi,i,maxi));
        colvec.append(index0fromnm(maxi,i,maxi));
        datvec.append(-bdr.gn(maxi,i,params));
        rowvec.append(index0fromnm(maxi,i,maxi));
        colvec.append(index0fromnm(maxi,i+1,maxi))#+(maxi+1));
        datvec.append(bdr.dm(maxi,i+1,params));
        rowvec.append(index0fromnm(maxi,j,maxi));
        colvec.append(index0fromnm(maxi,j,maxi)-1);
        datvec.append(bdr.bn(maxi-1,j,params));
        #rowvec.append(index0fromnm(maxi,i,maxi));
        #colvec.append(index0fromnm(maxi,i,maxi)-(maxi+1)+1);
        #datvec.append(bdr.bm(maxi+1,i-1,params));
    #full full
    rowvec.append(index0fromnm(maxi,maxi,maxi));
    colvec.append(index0fromnm(maxi,maxi,maxi));
    datvec.append(-bdr.gnm(maxi,maxi,params));
    rowvec.append(index0fromnm(maxi,maxi,maxi));
    colvec.append(index0fromnm(maxi,maxi,maxi)-1);
    datvec.append(bdr.bn(maxi-1,maxi,params));
    trix = sparse.coo_matrix((datvec,(rowvec,colvec)), shape=((maxi+1)*(maxi+1),(maxi+1)*(maxi+1)))
    return trix

def sparsepinit(consts,maxi):
    S0=int(ceil(consts[3]*consts[1]/consts[2]));
    I0=int(ceil(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]));
    initcond=((I0 - 1)*maxi + S0-1)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def sparseJempinit(consts,maxi):
    S0=int(ceil(consts[3]*consts[1]/consts[2]));
    I0=int(ceil(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]));
    initcond=((I0 - 1)*(maxi) + S0)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=((maxi)*(maxi+1),1));
    return pinit

def sparseJemP(consts,maxi,S0,I0):
    #S0=int(ceil(consts[3]*consts[1]/consts[2]));
    #I0=int(ceil(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]));
    initcond=((I0 - 1)*(maxi) + S0)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=((maxi)*(maxi+1),1));
    return pinit

def sparsevec(xIC,yIC,maxi):
    initcond=((yIC - 1)*maxi + xIC-1)
    rowvec = [initcond]; colvec = [0]; datvec = [1.0];
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def exitvec(a,kc,maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi):
                rowvec.append(i);
                colvec.append(0);
                datvec.append(-bdr.dm(1,i+1,a,kc));
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

def negones(maxi):
    rowvec = []; colvec = []; datvec = [];
    for i in range(maxi*maxi):
                rowvec.append(i);
                colvec.append(0);
                datvec.append(-1.0);
    pinit = sparse.coo_matrix((datvec,(rowvec,colvec)),shape=(maxi*maxi,1));
    return pinit

