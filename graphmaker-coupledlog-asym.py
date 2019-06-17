import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
from scipy.optimize import curve_fit
import numpy as np

koverk=5;#really this is 10*K1/K2 (or K2/K1)

#listofKs=[j for j in range(8,128+8,8)];
listofKs=[j for j in range(8,56+4,4)];
#listofas=[0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
listofas=[0.,.1,.2,.3,.4,.5]#listofas=[0.,.2,.4,.6,.8,1.]#listofas=[tena/10. for tena in range(1,koverk)];
numberofas=len(listofas);numberofKs=len(listofKs);

tdatumK = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
#tdatumtime = [[-1. for j in range(numberofKs)] for i in range(numberofas)];
# file's structure is K, MTE

''' for regular K (4 to 64) or (2 to 256) '''
for hay in range(numberofas):
#  fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\brokensymmetry\\kramers-nosym-20180223-K2overK1is%ia12is10a21is%i.txt'%(koverk,int(10*listofas[hay])),"r")
  fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\brokensymmetry\\kramers-nosym-20180810-K2overK1is5a12overa21is4a21is%i.txt'%(int(10*listofas[hay])),"r")
#  fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\brokensymmetry\\kramers-nosym-20181212-K2overK1is1a12is5a21is%i.txt'%(int(10*listofas[hay])),"r")
  crows=0;
  #for line in fly1:
  for crows in range(numberofKs):
    #line = fly1[crows];
    line = fly1.readline();
    prt = line.split("\t");
    if(float(prt[1])>0 and np.log10(abs(float(prt[1])))<13.):
      tdatumK[hay].append(float(prt[0]));
      tdatumtime[hay].append(float(prt[1]));
      #tdatumtime[hay][crows]=float(prt[1]);
    crows=crows+1;
    #if(crows>=numberofKs):
    #    break;
#for line in fly1:
#    prt=line.split("\t");
#    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.):#1E13 is a reasonable cutoff beyond which the algo doesn't work
#        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
#        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[2]));
  fly1.close()
  del fly1

def logansatz(kay,eff,gee,ach):
#    dummy = ach + gee*np.log(kay) + eff*kay; #assumes K is a single value
    dummy = ach + [gee*np.log(k) for k in kay] + [eff*k for k in kay]; #assumes K is a list
    return dummy

floist = ['g' for i in range(numberofas)]
glist = ['g' for i in range(numberofas)]
hlist = ['g' for i in range(numberofas)]
floisterr = ['g' for i in range(numberofas)]; glisterr = ['g' for i in range(numberofas)]; hlisterr = ['g' for i in range(numberofas)]
for i in range(numberofas):
    #print(i,tdatumK[i], np.log(tdatumtime[i]))
    #popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]), bounds=([-100.,-2.,-5.], [+100.,+2.,+5.]))
    popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]))
    [floist[i],glist[i],hlist[i]] = popt;
    [floisterr[i],glisterr[i],hlisterr[i]] = np.sqrt(np.diag(pcov))
#print(listofas,floist)

fntsz=13;
"""
''' this is the symmetric graph, for comparison '''
limitind=[np.exp(k)/k/2. for k in tdatumK[-1]]
limitWFM=[np.log(2)*k for k in tdatumK[-1]]
''' now plot mean time conditioned on successful invasion - ummm maybe mistitled'''
tauvK=plt.figure(2)
plot=tauvK.add_subplot(111)###
#plt.semilogy(tdatumK[-1],limitind,'k--',label='indep. limit')
for i in range(numberofas):
  if i%2==0:
    plt.semilogy(tdatumK[i],tdatumtime[i],'-',label='$a$ = %.1f'%(i/float(numberofas-1)),linewidth=1)
    plt.semilogy(tdatumK[i],np.exp(logansatz(tdatumK[i],floist[i],glist[i],hlist[i])),'--',label='$a$ = %.1f fit'%(i/float(numberofas-1)),linewidth=1)
#plt.semilogy(tdatumK[-1],limitind,'k--',linewidth=2)#to layer it atop
#plt.semilogy(tdatumK[-1],limitWFM,'r:',label='Moran limit',linewidth=3)
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('mean fixation time',fontsize=fntsz+4)
plt.xlabel(r'carrying capacity $K$',fontsize=fntsz+4)
plt.xlim(0,132)
plt.ylim(1E0,1E8)#12)
#plt.legend(loc=2,fontsize=fntsz)#1 - 4 are the corners
plt.show()
"""
"""
''' now plot mean time conditioned on successful invasion vs $a$ - ummm maybe mistitled'''
plt.figure(5)
for i in range(len(listofKs)-1,0-1,-1):
  if i%4==0:
    plt.semilogy(adatumK[i],adatumtime[i],'-',label='$K$ = %s'%(listofKs[i]),linewidth=3)
plt.ylabel('mean invasion time',fontsize=fntsz)
plt.xlabel('niche overlap, $a$',fontsize=fntsz)
plt.xlim(0,1)
plt.ylim(1E0,1E8)
plt.legend(loc=2,fontsize=fntsz)#1 - 4 are the corners, 10 is middle middle
plt.show()
"""

''' now plot f,g vs $a$'''
#from Gaussian, [[0.,.2,.4,.6,.8,1.],[0.118421, 0.132854, 0.15303, 0.111471, 0.0313293, 0.]] for K1=K2, a1=.5
ansatzplot=plt.figure(7)
plot=ansatzplot.add_subplot(111)###
plt.plot([0,.5], [1,0], 'o', color='blue')#expected f(a)
plt.plot([0,.5], [-1,1], 'o', color='green')#expected g(a)
plt.plot([0,.5], [-np.log(1.),np.log(np.log(2.))], 'o', color='orange')#expected h(a)
#plt.plot([0.,.2,.4,.6,.8,1.],[0.118421, 0.132854, 0.15303, 0.111471, 0.0313293, 0.],'o',color='red',label='Gaussian')
plt.errorbar(listofas,floist,fmt='--',label='$f(a)$',yerr=floisterr, color='blue')
plt.errorbar(listofas,glist,fmt=':',label='$g(a)$',yerr=glisterr, color='green')
plt.errorbar(listofas,hlist,fmt='-.',label='$h(a)$',yerr=hlisterr, color='orange')
#plt.plot(listofas,flist,label='$f(a)$', color='blue')
#plt.plot(listofas,glist,label='$g(a)$', color='green')
#plt.plot(listofas,hlist,label='$h(a)$', color='orange')
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('parameter values',fontsize=fntsz+4)
plt.xlabel(r'niche overlap $a_{21}$',fontsize=fntsz+4)# ($a_{12}=0.5$)
#plt.xlim(0,1)
#plt.ylim(-2,2)
plt.legend(fontsize=fntsz+3)
plt.show()
#ansatzplot.savefig("asym-K1overK2is2a12overa21is4.pdf",bbox_inches='tight')
print(floist)
