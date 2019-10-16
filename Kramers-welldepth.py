# want to get the residence times of timemap, use them as a Kramers' potential, find the well depth, fit it vs K for different a's
# 24/07/2018

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


listofKs=range(16,132+1,4);listofas=[0.0,0.2,0.4,0.6,0.8,1.0];
#listofKs=[8,16,32,64,128]; listofas=[tena/10. for tena in range(0,10+1)];
numberofKs=len(listofKs); numberofas=len(listofas);

tdatum = [[0. for j in range(numberofKs)] for i in range(numberofas)];
#tdatum = [[0.,0.,0.,0.,0.] for i in range(numberofas)];

# file's structure is a matrix of res time for species1pop,species2pop
#for tena in range(0,10+1):
for tena in range(numberofas):
    for k in range(numberofKs):
        kay=listofKs[k]
        startint=int(kay*1./(1.+listofas[tena])); endint=kay;
        fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixexpmodularized\\timemaps\\timemap-a%i-K%i-20180730.txt' %(2*tena,kay),"r")
        counterspell=1; tempup=0;tempdown=0;
        for line in fly1:
            prt=line.split("\t");
            if(counterspell==startint):
                tempdown=float(prt[counterspell-1]);
            if(counterspell==endint):
                tempup=float(prt[0]);
                break #you have all you need from this file (starting<endint)
            counterspell=counterspell+1
        fly1.close()
#        print(tena,kay,tempup-tempdown)
        tdatum[tena][k]=abs(tempup-tempdown); #all but two of these values are negative?!?!
#print(tdatum)

def logansatz(kay,eff,gee,ach):
#    dummy = ach + gee*np.log(kay) + eff*kay; #assumes K is a single value
    dummy = ach + [gee*np.log(k) for k in kay] + [eff*k for k in kay]; #assumes K is a list
    return dummy

flist = ['g' for i in range(numberofas)]
glist = ['g' for i in range(numberofas)]
hlist = ['g' for i in range(numberofas)]
flisterr = ['g' for i in range(numberofas)]; glisterr = ['g' for i in range(numberofas)]; hlisterr = ['g' for i in range(numberofas)]
for i in range(numberofas):
    #popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]), bounds=([-100.,-2.,-5.], [+100.,+2.,+5.]))
    popt, pcov = curve_fit(logansatz, listofKs, np.log(tdatum[i]))
    [flist[i],glist[i],hlist[i]] = popt;
    [flisterr[i],glisterr[i],hlisterr[i]] = np.sqrt(np.diag(pcov))
print(listofas,flist)

fntsz=13;
''' now plot f,g (of well depth) vs $a$'''
hopethisworks=plt.figure(0)
plot=hopethisworks.add_subplot(111)###
plt.errorbar(listofas,flist,fmt='--',label='welldepth $f(a)$',yerr=flisterr, color='blue')
plt.errorbar(listofas,glist,fmt=':',label='welldepth $g(a)$',yerr=glisterr, color='green')
#plt.errorbar(listofas,hlist,fmt='-.',label='$h(a)$',yerr=hlisterr, color='orange')
#plt.plot(listofas,flist,label='$f(a)$', color='blue')
#plt.plot(listofas,glist,label='$g(a)$', color='green')
#plt.plot(listofas,hlist,label='$h(a)$', color='orange')
plt.plot([0,1], [1,0], 'o', color='blue')#expected f(a)
plt.plot([0,1], [-1,1], 'o', color='green')#expected g(a)
#plt.plot([0,1], [-np.log(2.),np.log(np.log(2.))], 'o', color='orange')#expected h(a)
plot.tick_params(axis='both', which='major', labelsize=fntsz)###
plot.tick_params(axis='both', which='minor', labelsize=fntsz)###
plt.ylabel('parameter values',fontsize=fntsz+4)
plt.xlabel(r'niche overlap $a$',fontsize=fntsz+4)
plt.xlim(-0.1,1.1)
plt.ylim(-1.5,1.4)
#plt.legend(fontsize=fntsz+3)
#plt.show()
#ansatzplot.savefig("functionalKa9-welldepth.pdf",bbox_inches='tight')

fly10 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-fudge6-all-king.txt',"r")
listofKs=[j for j in range(4,132+4,4)];listofas=[tena/10. for tena in range(0,10+1)];
numberofas=len(listofas);numberofKs=len(listofKs);
tdatumK = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
for line in fly10:
    prt=line.split("\t");
    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.-6.):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[2])*1E6);#this is what fudge6 means
def logansatz0(kay,eff,gee,ach):
    dummy = ach + [gee*np.log(k) for k in kay] + [eff*k for k in kay]; #assumes K is a list
    return dummy
flist0 = ['g' for i in range(numberofas)]
glist = ['g' for i in range(numberofas)]
hlist = ['g' for i in range(numberofas)]
flisterr = ['g' for i in range(numberofas)]; glisterr = ['g' for i in range(numberofas)]; hlisterr = ['g' for i in range(numberofas)]
for i in range(numberofas):
    popt, pcov = curve_fit(logansatz0, tdatumK[i], np.log(tdatumtime[i]))
    [flist0[i],glist[i],hlist[i]] = popt;
    [flisterr[i],glisterr[i],hlisterr[i]] = np.sqrt(np.diag(pcov))
plt.errorbar(listofas,flist0,fmt='--',label='$f(a)$',yerr=flisterr, color='yellow')
plt.errorbar(listofas,glist,fmt=':',label='$g(a)$',yerr=glisterr, color='orange')
#plt.plot(listofas,flist0,label='$f(a)$', color='yellow')
#plt.plot(listofas,glist,label='$g(a)$', color='green')
#plt.plot(listofas,hlist,label='$h(a)$', color='orange')
plt.legend(fontsize=fntsz+3)
plt.show()
hopethisworks.savefig("welldepth.pdf",bbox_inches='tight');
fly10.close()
