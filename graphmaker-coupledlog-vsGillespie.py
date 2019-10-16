import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.cm as cm

fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-fudge6-all-king.txt',"r")
#fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-20180624.txt',"r")

listofKs=[j for j in range(4,132+4,4)];listofas=[tena/10. for tena in range(0,10+1)];
#listofKs=[j for j in range(16,64+4,8)];listofas=[tena/20. for tena in range(-2,24+1)];
numberofas=len(listofas);numberofKs=len(listofKs);

tdatumK = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
adatumK = [[] for i in range(numberofKs)];
adatumtime = [[] for i in range(numberofKs)];
# file's structure is K, a, MTE

''' for regular K (4 to 64) or (2 to 256) '''
for line in fly1:
    prt=line.split("\t");
#    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.-6.):
    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.):#1E13 is a reasonable cutoff beyond which the algo doesn't work
#    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[0]) not in tdatumK[int((numberofas-1)*float(prt[1]))]):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[2])*1E6);#this is what fudge6 means
#        tdatumK[int(10*float(prt[1])*2+2)].append(float(prt[0]));
#        tdatumtime[int(10*float(prt[1])*2+2)].append(float(prt[2]));#this is what fudge6 means
        adatumK[listofKs.index(int(prt[0]))].append(float(prt[1]));
        adatumtime[listofKs.index(int(prt[0]))].append(float(prt[2])*1E6);#this is what fudge6 means
fly1.close()

tGillK = [[] for i in range(numberofas)];
tGilltime = [[] for i in range(numberofas)];
fly2 = open('C:\\Users\\lenov\\Documents\\python\\tvKGillespie-20180808.txt',"r")
for line in fly2:
    prt=line.split("\t");
    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.):#1E13 is a reasonable cutoff beyond which the algo doesn't work
        tGillK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tGilltime[int((numberofas-1)*float(prt[1]))].append(float(prt[2]));
fly2.close()

fntsz=13;
cm_subsection = np.linspace(0.0,0.8,numberofas)
colors = [ cm.viridis(x) for x in cm_subsection ]

limitind=[np.exp(k)/k/2. for k in tdatumK[-1]]
limitWFM=[np.log(2)*k for k in tdatumK[-1]]
''' now plot mean time conditioned on successful invasion - ummm maybe mistitled'''
tauvK=plt.figure(2)
plot=tauvK.add_subplot(111)###
#plt.semilogy(tdatumK[-1],limitind,'k--',label='indep. limit')
for i in range(numberofas):
  if i%2==0:
    #plt.semilogy(tdatumK[i],tdatumtime[i],'-',color=cm.brg(listofas[i]),label='$a$ = %.1f'%(i/float(numberofas-1)))
    plt.semilogy(tdatumK[i],tdatumtime[i],marker='.',ls=":",c=colors[i],label='$a$ = %.1f'%(i/float(numberofas-1)))
    plt.semilogy(tGillK[i],tGilltime[i],marker='+',ls="--",c=colors[i])#,label='$a$ = %.1f'%(i/float(numberofas-1)),linewidth=3)
#plt.semilogy(tdatumK[-1],limitind,'k--',linewidth=2)#to layer it atop
#plt.semilogy(tdatumK[-1],limitWFM,'r:',label='Moran limit',linewidth=3)
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('mean fixation time',fontsize=fntsz+4)
plt.xlabel(r'carrying capacity $K$',fontsize=fntsz+4)
plt.xlim(4,20)
plt.ylim(1E0,1E5)
plt.legend(loc=1,fontsize=fntsz)#1 - 4 are the corners
plt.show()
#tauvK.savefig("coupled-logistic-data-vs-Gillespie.pdf",bbox_inches='tight')
