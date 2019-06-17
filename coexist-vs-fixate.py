import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
from scipy.optimize import curve_fit
import numpy as np

fly1 = open('C:\\Users\\lenov\\Documents\\python\\matrixsparse2\\kramers-tauvKa-fudge6-all-king.txt',"r")

listofKs=[j for j in range(4,132+4,4)];listofas=[tena/10. for tena in range(0,10+1)];
numberofas=len(listofas);numberofKs=len(listofKs);

tdatumK = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
adatumK = [[] for i in range(numberofKs)];
adatumtime = [[] for i in range(numberofKs)];
# file's structure is K, a, MTE

''' for K from 4 to 132 '''
for line in fly1:
    prt=line.split("\t");
    if(float(prt[2])>0 and np.log10(abs(float(prt[2])))<13.-6.):
#    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[0]) not in tdatumK[int((numberofas-1)*float(prt[1]))]):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[2])*1E6);#this is what fudge6 means
        adatumK[listofKs.index(int(prt[0]))].append(float(prt[1]));
        adatumtime[listofKs.index(int(prt[0]))].append(float(prt[2])*1E6);#this is what fudge6 means
fly1.close()

def logansatz(kay,eff,gee,ach):
#    dummy = ach + gee*np.log(kay) + eff*kay; #assumes K is a single value
    dummy = ach + [gee*np.log(k) for k in kay] + [eff*k for k in kay]; #assumes K is a list
    return dummy

flist = ['g' for i in range(numberofas)]
glist = ['g' for i in range(numberofas)]
hlist = ['g' for i in range(numberofas)]
flisterr = ['g' for i in range(numberofas)]; glisterr = ['g' for i in range(numberofas)]; hlisterr = ['g' for i in range(numberofas)]
for i in range(numberofas):
    popt, pcov = curve_fit(logansatz, tdatumK[i], np.log(tdatumtime[i]))
    [flist[i],glist[i],hlist[i]] = popt;
    [flisterr[i],glisterr[i],hlisterr[i]] = np.sqrt(np.diag(pcov))

def coexistenceDifference(kay,aye):#returns the difference between Moran and the data; if positive then fixation, else coexistence
    moranTime=np.log(np.log(2)*kay);
    trueTime=np.interp(aye,listofas,hlist)+np.log(kay)*np.interp(aye,listofas,glist)+kay*np.interp(aye,listofas,flist);
    return moranTime-trueTime;
def coexistenceSign(kay,aye):#returns the sign of the difference of Moran and the data; if positive then fixation, else coexistence
    moranTime=np.log(np.log(2)*kay);
    trueTime=np.interp(aye,listofas,hlist)+np.log(kay)*np.interp(aye,listofas,glist)+kay*np.interp(aye,listofas,flist);
    return np.sign(moranTime-trueTime);
def coexistenceOtherSign(kay,aye):#returns the sign of the difference of the data and Moran; if NEGATIVE then fixation, else coexistence
    moranTime=np.log(np.log(2)*kay);
    trueTime=np.interp(aye,listofas,hlist)+np.log(kay)*np.interp(aye,listofas,glist)+kay*np.interp(aye,listofas,flist);
    return (np.sign(trueTime-moranTime)+1.)/4.;
def doesItCoexist(kay,aye):#returns a boolean, whether Moran time is longer (ie fixates) or real time is longer (ie coexists)
    moranTime=np.log(np.log(2)*kay);
    trueTime=np.interp(aye,listofas,hlist)+np.log(kay)*np.interp(aye,listofas,glist)+kay*np.interp(aye,listofas,flist);
    if(moranTime>trueTime):
        return False
    elif(moranTime<trueTime):
        return True
    else:
        print("something has gone terribly wrong")
def toANumber(kay,aye):
    if(doesItCoexist(kay,aye)):
        return 1;
    else:
        return 0;

fntsz=13;
fntsz=17; #inset!

num = 1000; delta = 1./num;
#aaa = np.arange(-0.7, 1.0+.000001, delta)
#kkk = np.arange(48.,1.-0.000001, -delta)#this has to be reversed to imshow correctly - cf. contour plot which is right
aaa = np.arange(1.0,0.8-.000001, -delta)#flipit
kkk = np.arange(1.,25.+0.000001, +delta)#flipit
#X, Y = np.meshgrid(aaa,kkk)#unflip
Y, X = np.meshgrid(kkk,aaa)#flipit
#Z = toANumber(Y, X)
#Z = coexistenceDifference(Y, X)
Z = coexistenceOtherSign(Y, X)
"""
pudding=plt.figure(1)
plt.ylabel('carrying capacity $K$',fontsize=fntsz)
plt.xlabel('niche overlap $a$',fontsize=fntsz)
CS = plt.contourf(X, Y, Z)
plt.clabel(CS, inline=1, fontsize=fntsz)
plt.show()

firaxis=np.full(len(aaa),0.)
topaxis=np.full(len(aaa),150)
regionfig,(ax) = plt.subplots()
ax.fill_between(aaa,firaxis,topaxis,where=Z>0,interpolate=True)
plt.ylabel(r'niche overlap $a$',fontsize=fntsz+4)
plt.xlabel(r'carrying capacity $K$',fontsize=fntsz+4)
ax.set_aspect(aspect='equal')
ax.plot([0], [0], 'o', color='grey')
plt.xlim(0,1.)
plt.ylim(4,128)
plt.show()
"""

''' a and K; black is fixation '''
def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
figgypud=plt.figure(3)
ax = figgypud.add_subplot(111,aspect='equal')
#extent = np.min(aaa), np.max(aaa), np.min(kkk), np.max(kkk)#unflip
extent = np.min(kkk), np.max(kkk), np.min(aaa), np.max(aaa)#flipit
ax.imshow(Z,cmap=plt.cm.gray,extent=extent)
#im1 = plt.imshow(Z, cmap=plt.cm.gray, interpolation='nearest',extent=extent)
#plt.imshow(Z,cmap='hot',aspect='equal')
plt.box(on=None)
plt.ylim(1.0,0.801)
plt.xlim(1,25)
plt.xlabel('carrying capacity $K$',fontsize=fntsz)
plt.ylabel('niche overlap $a$    ',fontsize=fntsz)
forceAspect(ax,aspect=2.2)
plt.show()
#figgypud.savefig("coexist-vs-fixate27.jpeg",bbox_inches='tight',transparent=True)
