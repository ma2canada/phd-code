''' companion plotting file to kramers2-timemap.py '''

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.collections import LineCollection
import matplotlib.ticker as ticker
from scipy.integrate import solve_ivp

"""graphing deets"""
fig, ax = plt.subplots()
for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_locator(ticker.MaxNLocator(integer=True))

"""constants and reading the data"""
consts=[1,200,20,10];
NN=consts[1]; buffr=2; cutoff=NN*buffr; stopgraph=NN+5;
'''#this (below) is for tau vs init conds
fly1 = open('exttime-timemap-20200608.txt',"r")
justx = [i for i in range(70,150+1)]
justy = [i for i in range(1,20+1)]
plt.title("extinction time given initial conditions")
'''
fly1 = open('timemap-20200529.txt',"r")
justx = [i for i in range(cutoff+1)]
justy = [i+1 for i in range(cutoff)]
plt.title("residence times")
#'''#this (above) is for residence times
likemysoul=[];countr=0;
for line in fly1:
    likemysoul.append([]);
    prt=line.split("\t");
    for elem in prt:
        if(elem!='\n'):# or len(elem)!=0 or elem==prt[-1]):
            #print(elem,"okay?")
            likemysoul[-1].append(float(elem));
        else:
            print("this last element is just slash n")
    countr+=1;
likemysoul=np.array(likemysoul)
#print("argmax",np.argmax(np.array(likemysoul)))
print("shape",np.shape(np.array(likemysoul)))
print("x,y",len(justx),len(justy))
print("I hope it's x,y = ",len(likemysoul[0]),len(likemysoul))
print("each line was %i long, and there were %i lines" %(len(prt),countr))

"""plotting the data in a heatmap"""
plt.contourf(justx, justy, likemysoul, cmap=plt.cm.gray)#, norm=LogNorm())#,levels=25
#plt.contourf(likemysoul)
#plt.title("extinction time given initial conditions")
plt.xlabel("number susceptible (S)")
plt.ylabel("number infected (I)")
plt.xlim(75,135)
plt.ylim(0,20)
#plt.xlim(0,stopgraph)
#plt.ylim(0,stopgraph)
#plt.Axes.set_aspect(aspect='equal')
#plt.scatter(100,10,s=10,c='red')
#plt.show()

fly1.close()

"""trying to find a path in the heatmap"""
#path array
patharray=[]
#initcond
S0=int(np.ceil(consts[3]*consts[1]/consts[2]));
I0=int(np.ceil(consts[0]*(consts[2]-consts[3])*consts[1]/consts[2]/consts[3]));
patharray.append([S0,I0]);
hereS=(patharray[-1][0]); hereI=(patharray[-1][1]);
counter=0;
"""#this (below) is for tau vs init conds
while(hereS>=75 and hereI>0 and hereS<135 and hereI<20 and counter<1000):
    print(hereS,hereI)
    hereS=(patharray[-1][0]); hereI=(patharray[-1][1]);
    #default to left-up with restime
    if([hereS-1,hereI+1] not in patharray):
        nextS=hereS-1;nextI=hereI+1;
        resttime=likemysoul[nextI-1][nextS-75] #-1 to get to 0-indexing
    else:
        resttime=-1.0
#    print(resttime)
    #check left, right, down
    if(likemysoul[hereI-1][hereS-75-1]>resttime and [hereS-1,hereI] not in patharray):
        nextS=hereS-1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS-75-1]
    if(likemysoul[hereI-1][hereS-75+1]>resttime and [hereS+1,hereI] not in patharray):
        nextS=hereS+1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS-75+1]
    if(likemysoul[hereI-1-1][hereS-75]>resttime and [hereS,hereI-1] not in patharray):
        nextS=hereS;nextI=hereI-1;
        resttime=likemysoul[hereI-1-1][hereS-75]
    patharray.append([nextS,nextI]);
    counter+=1;
"""
while(hereS>=0 and hereI>0 and hereS<cutoff-1 and hereI<cutoff-1 and counter<163):
    print(hereS,hereI)
    hereS=(patharray[-1][0]); hereI=(patharray[-1][1]);
    #default to left-up with restime
    if([hereS-1,hereI+1] not in patharray):
        nextS=hereS-1;nextI=hereI+1;
        resttime=likemysoul[nextI-1][nextS] #-1 to get to 0-indexing
    else:
        resttime=-1.0
#    print(resttime)
    #check left, right, down
    if(likemysoul[hereI-1][hereS-1]>resttime and [hereS-1,hereI] not in patharray):
        nextS=hereS-1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS-1]
    if(likemysoul[hereI-1][hereS+1]>resttime and [hereS+1,hereI] not in patharray):
        nextS=hereS+1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS+1]
    if(likemysoul[hereI-1-1][hereS]>resttime and [hereS,hereI-1] not in patharray):
        nextS=hereS;nextI=hereI-1;
        resttime=likemysoul[hereI-1-1][hereS]
    patharray.append([nextS,nextI]);
    counter+=1;
#"""#this (above) is for residence times
print("number of steps was %i"%counter)
scats=(np.array(patharray)).T
points = np.array([scats[0], scats[1]]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
mynorm = plt.Normalize(0,counter)
lc = LineCollection(segments, cmap='spring', norm=mynorm)
lc.set_array(np.linspace(0, counter, counter+1))
ax.add_collection(lc)

"""trying to find a path in the heatmap in the reverse way"""
#path array
patharray=[]
#initcond
S0=100
I0=1
patharray.append([S0,I0]);
hereS=(patharray[-1][0]); hereI=(patharray[-1][1]);
counter=0;
while(hereS>=0 and hereI>0 and hereS<cutoff-1 and hereI<cutoff-1 and counter<36):
    print(hereS,hereI)
    hereS=(patharray[-1][0]); hereI=(patharray[-1][1]);
    #default to left-up with restime
    if([hereS-1,hereI+1] not in patharray):
        nextS=hereS-1;nextI=hereI+1;
        resttime=likemysoul[nextI-1][nextS] #-1 to get to 0-indexing
    else:
        resttime=-1.0
#    print(resttime)
    #check left, right, down
    if(likemysoul[hereI-1][hereS-1]>resttime and [hereS-1,hereI] not in patharray):
        nextS=hereS-1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS-1]
    if(likemysoul[hereI-1][hereS+1]>resttime and [hereS+1,hereI] not in patharray):
        nextS=hereS+1;nextI=hereI;
        resttime=likemysoul[hereI-1][hereS+1]
    if(likemysoul[hereI-1-1][hereS]>resttime and [hereS,hereI-1] not in patharray):
        nextS=hereS;nextI=hereI-1;
        resttime=likemysoul[hereI-1-1][hereS]
    patharray.append([nextS,nextI]);
    counter+=1;
print("number of steps was %i"%counter)
scats=(np.array(patharray)).T
points = np.array([scats[0], scats[1]]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
mynorm = plt.Normalize(0,counter)
lc = LineCollection(segments, cmap='winter', norm=mynorm)
lc.set_array(np.linspace(0, counter, counter+1))
ax.add_collection(lc)
#plt.plot(scats[0],scats[1],'r-')

eps=.1;
"""now to plot the time forward solution"""
stoptime=5.0
def SIR_model(t, y): return [consts[0]*consts[1] - consts[0]*y[0] - consts[2]*y[0]*y[1]/consts[1], consts[2]*y[0]*y[1]/consts[1] - consts[3]*y[1]]
sol = solve_ivp(SIR_model, [0.0, stoptime], [200-eps,eps],dense_output=True)
#print(sol.t)
#print(sol.y)
t = np.linspace(0, stoptime, 300)
x = sol.sol(t)
plt.plot(x[0],x[1],'r')
"""now to plot the time backward solution"""
stoptimeH=0.5
def SIR_Hamiltonian(t, y): return [consts[0]*(consts[1]-y[0])-consts[2]*y[0]*y[1]*y[3]/consts[1], consts[2]*y[0]*y[1]*(2*y[3]-y[2])/consts[1]-consts[3]*y[1], consts[0]*(y[2]-1)+consts[2]*(y[2]-y[3])*y[1]*y[3]/consts[1], consts[3]*(y[3]-1)+consts[2]*(y[2]-y[3])*y[0]*y[3]/consts[1]]
solH = solve_ivp(SIR_Hamiltonian, [0.0, stoptimeH], [100,10,1-eps,1],dense_output=True)
#solH = solve_ivp(SIR_Hamiltonian, [0.0, stoptimeH], [200-eps,1+eps,1,1],dense_output=True)
#print(solH.t)
#print(solH.y)
tH = np.linspace(0, stoptimeH, 300)
xH = solH.sol(tH)
plt.plot(xH[0],xH[1],'g--',linewidth=3)

plt.show()
#print(patharray)
