import matplotlib.pyplot as plt
#import matplotlib.ticker as plticker
from numpy import log,linspace,exp
from scipy.integrate import odeint
#import coupledlogistic-integrator as cli
import loginvasion as liar
from matplotlib import cm

fly2 = open('C:\\Users\\lenov\\Documents\\scinet\\invasion maybe\\tauvKa-fiftyfifty-20170209-lowK.txt',"r")
#fly = open('C:\\Users\\lenov\\Documents\\scinet\\invasion maybe\\tauvKa-fiftyfifty-20161230-lin.txt',"r")
fly = open('C:\\Users\\lenov\\Documents\\scinet\\invasion maybe\\tauvKa-20161230-exp-long.txt',"r")
#fly = open('C:\\Users\\lenov\\Documents\\scinet\\invasion maybe\\tauvKa-fiftyfifty-20161230-lin-long.txt',"r")
fly3 = open('C:\\Users\\lenov\\Documents\\scinet\\invasion maybe\\tauvKa-fiftyfifty-20170210-higha.txt',"r")

numberofas=11;
listofKs=[str(2**(j+1)) for j in range(8)];

tdatumK = [[] for i in range(numberofas)];
tdatumprob = [[] for i in range(numberofas)];
tdatumtime = [[] for i in range(numberofas)];
tdatumtime2 = [[] for i in range(numberofas)];
adatumK = [[] for i in range(len(listofKs))];#I realize I could just take the previous arrays and transpose them; what I do here is inefficient
adatumprob = [[] for i in range(len(listofKs))];
adatumtime = [[] for i in range(len(listofKs))];
adatumtime2 = [[] for i in range(len(listofKs))];
# file's structure is K, a, n, m, invasion prob, failure prob, invasion time, fail time

''' and for low K  (2 to 16)'''
# file's structure is K, a, n, m, invasion prob, failure prob, invasion time, fail time
for line in fly2:
    prt=line.split("\t");
    if(prt[2]==prt[0] and prt[3]=='1'):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumprob[int((numberofas-1)*float(prt[1]))].append(float(prt[4]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[6]));
        tdatumtime2[int((numberofas-1)*float(prt[1]))].append(float(prt[7]));

''' for regular K (4 to 64) or (2 to 256) '''
for line in fly:
    prt=line.split("\t");
#    if(prt[2]==prt[0] and prt[3]=='1'):
    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[0]) not in tdatumK[int((numberofas-1)*float(prt[1]))]):
        tdatumK[int((numberofas-1)*float(prt[1]))].append(float(prt[0]));
        tdatumprob[int((numberofas-1)*float(prt[1]))].append(float(prt[4]));
        tdatumtime[int((numberofas-1)*float(prt[1]))].append(float(prt[6]));
        tdatumtime2[int((numberofas-1)*float(prt[1]))].append(float(prt[7]));
    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[1])<0.8):#.8 is where high $a$ starts
        adatumK[listofKs.index(prt[0])].append(float(prt[1]));
        adatumprob[listofKs.index(prt[0])].append(float(prt[4]));
        adatumtime[listofKs.index(prt[0])].append(float(prt[6]));
        adatumtime2[listofKs.index(prt[0])].append(float(prt[7]));

''' for high a '''
for line in fly3:
    prt=line.split("\t");
    if(prt[2]==prt[0] and prt[3]=='1' and float(prt[1]) not in adatumK[listofKs.index(prt[0])]):
        adatumK[listofKs.index(prt[0])].append(float(prt[1]));
        adatumprob[listofKs.index(prt[0])].append(float(prt[4]));
        adatumtime[listofKs.index(prt[0])].append(float(prt[6]));
        adatumtime2[listofKs.index(prt[0])].append(float(prt[7]));

def esseye(i,kay):#prod d/b
    Si=1.
    for j in range(i):
        Si=Si*(j+1)/kay
    return Si
def queeye(i,kay):# 1/d prod b/d
    qi=kay/i**2
    for j in range(i-1):
        qi=qi*kay/(j+1)
    return qi
def indepinvasen(n,kay):
    #probability of failure/extinction
    numer=0; denompart=0;
    for j in range(n,int(kay/2)+1):
        numer+=esseye(j+1,kay)
    for j in range(int(kay/2)+1):
        denompart+=esseye(j+1,kay)
    return numer/(1+denompart)
def indepinvas(kay):
    #probability of failure/extinction
    numer=0;
    for j in range(int(kay/2)+1):
        numer+=esseye(j+1,kay)
    return numer/(1+numer)
#    return 1/(1+numer)

def Morantausucc(kay):
    return kay*(kay-1)*log(kay/(kay-1))
#    return kay/2*(kay/2-1)*log(kay/2/(kay/2-1))
def Morantaufail(kay):
    return (kay-2)/kay*(log(kay)-(kay-1)*log(kay/(kay-1)))
#    return (kay/2-2)/kay/2*(log(kay/2)-(kay/2-1)*log(kay/2/(kay/2-1)))

def phi1(kay):
    numer=0.;
    for j in range(1,int(kay/2)+1):
        for i in range(1,j+1):
            numer+=indepinvasen(i,kay)*queeye(i,kay)
    denompart=0;
    for j in range(1,int(kay/2)+1):
        denompart+=esseye(j,kay)
    return numer/(1+denompart)
def antiphi1(kay):
    numer=0.;
    for j in range(1,int(kay/2)+1):
        for i in range(1,j+1):
            numer+=(1-indepinvasen(i,kay))*queeye(i,kay)
    denompart=0;
    for j in range(1,int(kay/2)+1):
        denompart+=esseye(j,kay)
    return numer/(1+denompart)
def indeptaufail(kay):#untouched but it looks like it's written correct but the plot is bad
    numer=0.;
    for j in range(int(kay)+1-1):
        for i in range(1,j+1+1):
            numer+=indepinvasen(i,kay)*queeye(i,kay)
    denompart=0;
    for j in range(int(kay)+1-1):
        denompart+=esseye(j+1,kay)
    return numer/indepinvas(kay)/(1+denompart)
def indeptausucc(kay):
    return antiphi1(kay)/(1-indepinvas(kay));

def smallKextinction(kay,aay):
    return kay/(1+aay*(kay-1))
def smallKextinction2(kay,aay):
    return 1/((1+aay*(kay-1))/kay + 1.)

t=linspace(0.,100.,1000)
def xvecdot(xvec,t,kay,aay):
    r1=1.0;K1=kay;a1=aay;r2=r1;K2=K1;a2=a1;
    return [r1*xvec[0]*(1-(xvec[0]+a1*xvec[1])/K1),
            r2*xvec[1]*(1-(xvec[1]+a2*xvec[0])/K2)]
def taufinder(K1,a1):
    [x0,y0] = [K1-1,0+1]
    solution=odeint(xvecdot,[x0,y0],t,args=(K1,a1))
    indx=0; currenty=y0;
    while currenty < K1/(1+a1)-1:
        indx+=1
        currenty = solution[indx,1]
    return t[indx]
Klist=[2*i for i in range(1,64)]
Tlist=[]
Tlist2=[]
for k in Klist:
    K1=k
    Tlist.append(taufinder(k,0.0))
    Tlist2.append(taufinder(k,0.5))

def kimura_s(kay,aay):
    return (kay-1+aay)/(aay*kay+1-aay)
def kimura_prob(s,N):
    return (1.-exp(-2*s))/(1.-exp(-4*N*s))

wordsize=16;
number_of_lines= numberofas
#cm_subsection = linspace(0.0,1.0,number_of_lines)
#colors = [ cm.winter(x) for x in cm_subsection ]
cm_subsection = linspace(0.0,0.8,number_of_lines)
colors = [ cm.viridis(x) for x in cm_subsection ]

''' now plot probability of invasion '''
probvK=plt.figure(1)
plot=probvK.add_subplot(111)###
plt.plot([i for i in range(3,250)],[1-indepinvas(i) for i in range(3,250)],color='purple',ls='-',label='independent limit')
#plt.plot(tdatumK[0],[indepinvas(tdatumK[0][j]) for j in range(len(tdatumK[0]))],'k.',label='independent limit')
#plt.plot(tdatumK[-1],[liar.probvec2(tdatumK[-1][j])[0] for j in range(len(tdatumK[-1]))],'m--',label='trytry')#tdatumK[-1],
for i in range(numberofas):
  if i%2==0:
    tempa=i/float(numberofas-1);
    plt.plot(tdatumK[i],tdatumprob[i],marker='.',ls=":",label='$a$ = %.1f'%(i/float(numberofas-1)),color=colors[i])
#    plt.plot(tdatumK[i],[kimura_prob(kimura_s(tdatumK[i][j],tempa),tdatumK[i][j]) for j in range(len(tdatumK[i]))],'r--',label='Kimura')
plt.plot(tdatumK[-1],[2/tdatumK[-1][j] for j in range(len(tdatumK[-1]))],'g-',label='Moran limit')
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('invasion probability',fontsize=wordsize)
plt.xlabel('carrying capacity $K$',fontsize=wordsize)
plt.xlim(0,100)
plt.ylim(0,1)
plt.legend()
plt.show()
probvK.savefig("fiftyfifty-probvK.pdf",bbox_inches='tight')

''' now plot mean time conditioned on successful invasion '''
invvK=plt.figure(2)
plot=invvK.add_subplot(111)###
plt.plot(tdatumK[-1],[Morantausucc(tdatumK[-1][j]) for j in range(len(tdatumK[-1]))],'g-',label='Moran limit')
for i in range(numberofas-1,0-1,-1):
  if i%2==0:
    plt.plot(tdatumK[i],tdatumtime[i],marker='.',ls=":",label='$a$ = %.1f'%(i/float(numberofas-1)),color=colors[i])
#    if i!=0: plt.plot(tdatumK[i],[(mK**(i/float(numberofas-1))-1)/(i/float(numberofas-1)) for mK in tdatumK[i]],'-',label='(K^a-1)/a for $a$ = %.1f'%(i/float(numberofas-1)))
#plt.plot(tdatumK[-1],[indeptausucc(tdatumK[-1][j]) for j in range(len(tdatumK[-1]))],'k-',label='independent limit succ')
#plt.plot(tdatumK[-1],[indeptaufail(tdatumK[-1][j]) for j in range(len(tdatumK[-1]))],'k-',label='independent limit fail')
plt.plot(tdatumK[-1],[liar.tauinv(tdatumK[-1][j])[0] for j in range(len(tdatumK[-1]))],color='purple',ls='-',label='independent limit')
#plt.plot(Klist,Tlist,'-',label='$a=0.0$ deterministic')
#plt.plot(Klist,Tlist2,'-',label='$a=0.5$ deterministic')
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('mean invasion time',fontsize=wordsize)
plt.xlabel('carrying capacity $K$',fontsize=wordsize)
plt.xlim(0,100)
plt.ylim(0,40)
#plt.yscale('log')#!!!
#plt.xscale('log')
plt.legend(loc='upper left')
plt.show()
invvK.savefig("fiftyfifty-invtimevK.pdf",bbox_inches='tight')

''' now plot mean time conditioned on failure to invade '''
extvK=plt.figure(3)
plot=extvK.add_subplot(111)###
plt.plot(tdatumK[-1],[Morantaufail(tdatumK[-1][j]) for j in range(len(tdatumK[-1]))],'g-',label='Moran limit')
for i in range(numberofas-1,0-1,-1):
  if i%2==0:
    plt.plot(tdatumK[i],tdatumtime2[i],marker='.',ls=":",label='$a$ = %.1f'%(i/float(numberofas-1)),color=colors[i])
#plt.plot(tdatumK[-1],[indeptaufail(tdatumK[-1][j]) for j in range(len(tdatumK[-1]))],'k-',label='independent limit')
plt.plot(tdatumK[-1],[liar.taufai(tdatumK[-1][j])[0] for j in range(len(tdatumK[-1]))],color='purple',ls='-',label='independent limit')
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('mean failure time',fontsize=wordsize)
plt.xlabel('carrying capacity $K$',fontsize=wordsize)
plt.xlim(0,100)
plt.ylim(0,4)
#!plt.yscale('log')#!!!
plt.legend()
plt.show()
extvK.savefig("fiftyfifty-exttimevK.pdf",bbox_inches='tight')


''' now plot probability of invasion vs $a$'''
probva=plt.figure(4)
plot=probva.add_subplot(111)###
tempK=3;
#plt.plot([i/10 for i in range(0,11)],[1/(4./3.+i/10) for i in range(0,11)],'k.',label='small $K$')
#plt.plot([i/10 for i in range(0,11)],[tempK/(tempK+1+(tempK-1)*(i/10+tempK)+tempK-1+i/10) for i in range(0,11)],'k-.',label='alternate small $K$')
plt.plot(adatumK[-1],[1-adatumK[-1][j] for j in range(len(adatumK[-1]))],'m-',label='large $K$')
#for i in range(len(listofKs)-1,0-1,-1):
for i in [7,6,5,4,3,2,1,0]:
  if i%2==1:
    plt.plot(adatumK[i],adatumprob[i],marker='.',ls=":",label='$K$ = %s'%(listofKs[i]))
plt.plot([i/10 for i in range(0,10+1)],[tempK/(tempK+1+(tempK-1)*i/10) for i in range(0,10+1)],'k-',label='small $K$')
#plt.plot(adatumK[0],[1/(4./3.+adatumK[0][j]) for j in range(len(adatumK[0]))],'k.',label='small $K$')
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('invasion probability',fontsize=wordsize)
plt.xlabel('niche overlap $a$',fontsize=wordsize)
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend()
plt.show()
probva.savefig("fiftyfifty-probva.pdf",bbox_inches='tight')

''' now plot mean time conditioned on successful invasion vs $a$'''
invva=plt.figure(5)
plot=invva.add_subplot(111)###
#plt.plot(adatumK[-1],[smallKextinction(4,adatumK[-1][j]) for j in range(len(adatumK[-1]))],'k-.',label='$1/d_{mut}$')
#plt.plot([i/10 for i in range(0,11)],[1./(tempK/(1+(tempK-1)*i/10) + 1/1 + 1/(tempK-1) + tempK/(i/10+tempK-1)) for i in range(0,11)],'k-',label='stillwrong small $K$')
for i in range(len(listofKs)-1,0-1,-1):
  if i%2==1:
    plt.plot(adatumK[i],adatumtime[i],marker='.',ls=":",label='$K$ = %s'%(listofKs[i]))
plt.plot(adatumK[-1],[smallKextinction2(4,adatumK[-1][j]) for j in range(len(adatumK[-1]))],color='cyan',ls="-",label='$1/(d_{mut}+b_{mut})$')
#plt.plot(adatumK[-1],[smallKextinction(4,adatumK[-1][j]) for j in range(len(adatumK[-1]))],'r--',label='$1/(d_{mut})$')
#plt.plot(adatumK[1],adatumtime2[1],'-',label='$K$ = 4 failure') this is diff from the success time
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('mean invasion time',fontsize=wordsize)
plt.xlabel('niche overlap $a$',fontsize=wordsize)
plt.xlim(0,1)
#plt.ylim(0,4)
plt.legend()
plt.show()
invva.savefig("fiftyfifty-invtimeva.pdf",bbox_inches='tight')

''' now plot mean time conditioned on failure to invade vs $a$'''
extva=plt.figure(6)
plot=extva.add_subplot(111)###
#plt.plot(adatumK[-1],[smallKextinction(3,adatumK[-1][j]) for j in range(len(adatumK[-1]))],'k-.',label='$1/d_{mut}$')
#plt.plot([i/10 for i in range(0,11)],[1./(tempK/(1+(tempK-1)*i/10) + 1/1 + 1/(tempK-1) + tempK/(i/10+tempK-1)) for i in range(0,11)],'k-',label='stillwrong small $K$')
for i in range(len(listofKs)-1,0-1,-1):
  if i%2==1:
    plt.plot(adatumK[i],adatumtime2[i],marker='.',ls=":",label='$K$ = %s'%(listofKs[i]))
plt.plot(adatumK[-1],[smallKextinction2(4,adatumK[-1][j]) for j in range(len(adatumK[-1]))],c='cyan',ls="-",label='$1/(d_{mut}+b_{mut})$')
#plt.plot(adatumK[1],adatumtime[1],'-',label='$K$ = 4 success') this is different from the failure time
plot.tick_params(axis='both', which='major', labelsize=wordsize)###
plot.tick_params(axis='both', which='minor', labelsize=wordsize)###
plt.ylabel('mean failure time',fontsize=wordsize)
plt.xlabel('niche overlap $a$',fontsize=wordsize)
plt.xlim(0,1)
#plt.ylim(0,4)
plt.legend()
plt.show()
extva.savefig("fiftyfifty-exttimeva.pdf",bbox_inches='tight')

fly.close()
fly2.close()
fly3.close()
