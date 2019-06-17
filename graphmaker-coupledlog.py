import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import cm

fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-fudge6-all-king.txt',"r")#paper
#fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-20180624.txt',"r")

listofKs=[j for j in range(4,132+4,4)];listofas=[tena/10. for tena in range(0,10+1)];#paper
#listofKs=[j for j in range(16,64+4,8)];listofas=[tena/20. for tena in range(-2,24+1)];
numberofas=len(listofas);numberofKs=len(listofKs);

tdatumK = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
adatumK = [[] for i in range(numberofKs)];
adatumtime = [[] for i in range(numberofKs)];
# file's structure is K, a, MTE

"""
''' and for low K  (2 to 16)'''
# file's structure is K, a, n, m, invasion prob, failure prob, invasion time, fail time
for line in fly2:
    prt=line.split("\t");
    if(prt[2]==prt[0] and prt[3]=='1'):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumprob[int((numberofas-1)*float(prt[1]))].append(float(prt[4]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[6]));
        tdatumtime2[int((numberofas-1)*float(prt[1]))].append(float(prt[7]));
"""
''' for regular K (4 to 64) or (2 to 256) '''
for line in fly1:
    prt=line.split("\t");
    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.-6.):#paper
#    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.):#1E13 is a reasonable cutoff beyond which the algo doesn't work
#    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[0]) not in tdatumK[int((numberofas-1)*float(prt[1]))]):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));#paper
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[2])*1E6);#this is what fudge6 means#paper
#        tdatumK[int(10*float(prt[1])*2+2)].append(float(prt[0]));
#        tdatumtime[int(10*float(prt[1])*2+2)].append(float(prt[2]));#this is what fudge6 means
        adatumK[listofKs.index(int(prt[0]))].append(float(prt[1]));
        adatumtime[listofKs.index(int(prt[0]))].append(float(prt[2])*1E6);#this is what fudge6 means

def logansatz(kay,eff,gee,ach):
#    dummy = ach + gee*np.log(kay) + eff*kay; #assumes K is a single value
    dummy = ach + [gee*np.log(k) for k in kay] + [eff*k for k in kay]; #assumes K is a list
    return dummy

flist = ['g' for i in range(numberofas)]
glist = ['g' for i in range(numberofas)]
hlist = ['g' for i in range(numberofas)]
flisterr = ['g' for i in range(numberofas)]; glisterr = ['g' for i in range(numberofas)]; hlisterr = ['g' for i in range(numberofas)]
for i in range(numberofas):
    #print(i,tdatumK[i], np.log(tdatumtime[i]))
    #popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]), bounds=([-100.,-2.,-5.], [+100.,+2.,+5.]))
    popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]))
    [flist[i],glist[i],hlist[i]] = popt;
    [flisterr[i],glisterr[i],hlisterr[i]] = np.sqrt(np.diag(pcov))
print(listofas,flist)

fntsz=13;
cm_subsection = np.linspace(0.0,0.8,numberofas)
colors = [ cm.viridis(x) for x in cm_subsection ]

limitind=[np.exp(k)/k/2. for k in tdatumK[-1]]
limitWFM=[np.log(2)*k for k in tdatumK[-1]]
''' now plot mean time conditioned on successful invasion - ummm maybe mistitled'''
tauvK=plt.figure(2)
plot=tauvK.add_subplot(111)###
#plt.semilogy(tdatumK[-1],limitind,'k--',label='indep. limit')
plt.semilogy(tdatumK[-1],limitind,c='purple',ls='-',label='indep. limit')
for i in range(numberofas):
  if i%2==0:
    plt.semilogy(tdatumK[i],tdatumtime[i],marker='.',ls=":",c=colors[i],label='$a$ = %.1f'%(i/float(numberofas-1)))#,linewidth=3)
#plt.semilogy(tdatumK[-1],limitind,'k--',linewidth=2)#to layer it atop
plt.semilogy(tdatumK[-1],limitind,c='purple',ls='-')#to layer it atop
plt.semilogy(tdatumK[-1],limitWFM,'g-',label='Moran limit')#,linewidth=3)#'r:'
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('mean fixation time',fontsize=fntsz+4)
plt.xlabel(r'carrying capacity $K$',fontsize=fntsz+4)
plt.xlim(0,132)
plt.ylim(1E0,1E12)
plt.legend(loc=1,fontsize=fntsz)#1 - 4 are the corners
plt.show()
tauvK.savefig("coupled-logistic-data.pdf",bbox_inches='tight')

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
"""
''' now plot f,g vs $a$'''
ansatzplot=plt.figure(7)
plot=ansatzplot.add_subplot(111)###
plt.errorbar(listofas,flist,fmt='--',label='$f(a)$',yerr=flisterr, color='blue')
plt.errorbar(listofas,glist,fmt=':',label='$g(a)$',yerr=glisterr, color='green')
plt.errorbar(listofas,hlist,fmt='-.',label='$h(a)$',yerr=hlisterr, color='orange')
#plt.plot(listofas,flist,label='$f(a)$', color='blue')
#plt.plot(listofas,glist,label='$g(a)$', color='green')
#plt.plot(listofas,hlist,label='$h(a)$', color='orange')
plt.plot([0,1], [1,0], 'o', color='blue')#expected f(a)
plt.plot([0,1], [-1,1], 'o', color='green')#expected g(a)
plt.plot([0,1], [-np.log(2.),np.log(np.log(2.))], 'o', color='orange')#expected h(a)
#tempalist=np.asarray([0,0.2,0.4,0.6,0.8,1.0]);plt.plot(tempalist, (1.-tempalist)/2/(1.+tempalist), 'x', color='red')#FP Gaussian f(a)
tempalist=np.linspace(0,1,1000);plt.plot(tempalist, (1.-tempalist)/2/(1.+tempalist), color='red')#FP Gaussian f(a)
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('parameter values',fontsize=fntsz+4)
plt.xlabel(r'niche overlap $a$',fontsize=fntsz+4)
#plt.xlim(0,1)
#plt.ylim(-2,2)
plt.legend(fontsize=fntsz+3-5,loc=6)
plt.show()
#ansatzplot.savefig("functionalKa10.pdf",bbox_inches='tight')
"""
fly1.close()
