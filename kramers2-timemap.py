import maketrix as mk
from scipy.sparse.linalg import spsolve
from numpy import floor
import datetime
#K=32 takes ~3.5min [0r 7 for buffer 5?], K=64 takes ~1hour [40min?], 128 ~4 hrs

buffr = 2;
kconst = 64;
print(datetime.datetime.now());
for kconst in range(136,174,4):
  for aa in range(0,11,2):
    aconst = aa/10.
    homme = open("timemap-a%i-K%i-20180730.txt"%(aa,kconst),'w');
    n0=int(floor(kconst/(1.+aconst)));
    initcond=((n0 - 1)*kconst*buffr + n0-1)
    
    temp3=mk.makesparsetrix(aconst,kconst,kconst*buffr)
    print("matrix acquired at "+str(datetime.datetime.now()))
    pickyvec = mk.sparsepinit(aconst,kconst,kconst*buffr);
    sumthese = spsolve(temp3.tocsc(),pickyvec.tocsc());
    tau = -sum(sumthese);
    print("with kconst = %i (buffer %i) and aconst = %f you get tau = %f" %(kconst,buffr,aconst,tau)); print();
    for elem in range(len(sumthese)):
        homme.write("%f\t"%(-sumthese[elem]));
        if((elem+1)%(kconst*buffr)==0):
            homme.write("\n");
    homme.close()
print(datetime.datetime.now());

