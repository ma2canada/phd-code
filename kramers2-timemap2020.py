import maketrix2020 as mk20
from scipy.sparse.linalg import spsolve, eigs
import datetime
#N=200 takes ~3hrs, almost all of that time is making the matrix

buffr = 2;
consts = [1, 200, 20, 10]; #[mu, N, beta, Gamma]
maxtrix=consts[1]*buffr;
datemeplease=20200614;

print(datetime.datetime.now()); print();
###alotalot = open("exttime-timemap-%i.txt" %datemeplease,'w');
if(True):
### for ii in range(1,20+1):
###  for ss in range(70,150+1):
###    initcond=mk20.indexJfromnm(ss,ii,maxtrix)
###    """
    homme = open("timemap-%i.txt" %datemeplease,'w');
###    """
    
    temp3=mk20.makeJemtrix(consts,maxtrix)
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
    print("with the parameters [1, 200, 20, 10] you get tau = %f" %(tau)); print();
    #print("ya goof, the length is %i"%len(sumthese))
    hahaha=0
    for elem in range(len(sumthese)):
        hahaha=hahaha+1;
        homme.write("%f"%(-sumthese[elem]));
        if(mk20.nfromindexJ(elem,maxtrix)==maxtrix):
            homme.write("\n");
        else:
            homme.write("\t");
    homme.close()
    #print("each line contains %i elements, the last of which is %f" %(hahaha,-sumthese[-1]))
    filer = open("exttime-%i.txt" %datemeplease,'w');
    filer.write("%f"%(tau));
    filer.close()
    """
    alotalot.write("%f\t"%(tau));
    if(mk20.nfromindexJ(initcond,maxtrix)==150):
        alotalot.write("\n");
alotalot.close()
"""###
print(datetime.datetime.now());

'''for testing; archaic
temp4=mk20.maketrix_test(consts,maxtrix)
for line in temp4:
    print(line)
#print(temp4.toarray())
'''
