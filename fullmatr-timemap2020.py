import maketrix2020 as mk20
from scipy.sparse.linalg import spsolve, eigs, expm_multiply
from scipy.linalg import expm
from numpy import dot, where, sort, shape, array, flip
import datetime
#N=200 takes ~3hrs, almost all of that time is making the matrix

buffr = 2;
consts = [1, 200, 20, 10]; #[mu, N, beta, Gamma]
maxtrix=consts[1]*buffr;
datemeplease=20200624;

print(datetime.datetime.now()); print();
###alotalot = open("exttime-timemap-%i.txt" %datemeplease,'w');
if(True):
### for ii in range(1,20+1):
###  for ss in range(70,150+1):
###    initcond=mk20.indexJfromnm(ss,ii,maxtrix)
###    """
    homme = open("lasttraj-%i.txt" %datemeplease,'w');
###    """
    
    temp3=mk20.makeJemtrix(consts,maxtrix);
    #temp2=temp3.todense();
    print("matrix acquired at "+str(datetime.datetime.now()))#+" for initcond = ",initcond)
    print();
    '''
    print("smallest eigenvalues are:",eigs(temp3,which='SM')[0])
    print(" largest eigenvalues are:",eigs(temp3,which='LM')[0])
    print("the matrix itself is:")#,temp3.toarray())
    for line in temp3.toarray():
        print(line)
    '''
    pickyvec = mk20.sparseJempinit(consts,maxtrix);
###    pickyvec = mk20.sparseJemP(consts,maxtrix,mk20.nfromindexJ(initcond,maxtrix),mk20.mfromindexJ(initcond,maxtrix));
    sumthese = spsolve(temp3.tocsc(),pickyvec.tocsc());
    tau = -sum(sumthese);
###    """
    print("with the parameters [1, 200, 20, 10] you get tau = %f, at time "%(tau)+str(datetime.datetime.now())); print();
    #finalvec=dot(expm(tau*temp2),pickyvec.todense())
    finalvec=array(expm_multiply((tau*temp3),pickyvec).todense())
    #print(finalvec,"which is",type(finalvec),"and shape",shape(finalvec))
    bigstates=[];
    for elem in finalvec:
        bigstates.append(elem[0])
    print("final vector at "+str(datetime.datetime.now())); print();
    '''now to pick the most probable states'''
    bigstates=sort(array(bigstates))[-100:]
    #bigstates=sort(finalvec)[-100:]
    bigstates=flip(bigstates)
    #print(bigstates)
    hahaha=0
    for elem in bigstates:
        hahaha=hahaha+1;
        bigindex=where(finalvec==elem)[0][0]
        homme.write("%i\t%i" %(mk20.nfromindexJ(bigindex,maxtrix),mk20.mfromindexJ(bigindex,maxtrix)) );
        #print(mk20.nfromindexJ(bigindex,maxtrix),mk20.mfromindexJ(bigindex,maxtrix))
        if(hahaha!=len(bigstates)):
            homme.write("\n");
    homme.close()
    """
    alotalot.write("%f\t"%(tau));
    if(mk20.nfromindexJ(initcond,maxtrix)==150):
        alotalot.write("\n");
alotalot.close()
"""###
print(datetime.datetime.now());
