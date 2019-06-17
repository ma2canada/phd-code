'''
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Futura']
'''
import matplotlib.pyplot as plt
import numpy as np
#import pylab as pl
from matplotlib import rc
rc('text', usetex=True)

'''
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}
pl.rc('font', **font)

plt.style.use('seaborn-deep')
#font = {'family':'Freeserif','size':16, 'serif': ['computer modern roman']}
font = {'size':16, 'serif': ['computer modern roman']}
pl.rc('font',**font)
#plt.rc('font', serif=['Helvetica'])

from matplotlib.font_manager import FontProperties
font0 = FontProperties()
font0.set_family(serif=['Helvetica'])

import matplotlib as mpl
#font={"Helvetia"}
mpl.rc('font', serif=['Courrierasdf'])
'''

xdataneg = np.arange(-2,0,0.01)
xdatapos = np.arange(0,+1,0.01)
xdatatot = np.concatenate((xdataneg,xdatapos))
y1datneg = np.array([1. for i in range(len(xdataneg))])
y1datpos = np.array([1. for i in range(len(xdatapos))])
y1dattot = np.concatenate((y1datneg,y1datpos))
y2datneg = 1./xdataneg
y2datpos = np.array([-2 for i in range(len(xdatapos))])
y2dattot = np.concatenate((y2datneg,y2datpos))
xdatapos2 = np.arange(1,2,0.01)
y1dattot2 = np.array([1. for i in range(len(xdatapos2))])
'''
numpoints = 101
xdata = [-2.+(1--2)*i/(numpoints-1) for i in range(numpoints)]
y1dat = [1. for i in range(len(xdata))]
y2dat = [-1./xdata[i] for i in range(len(xdata))]
#y2dat = np.array([-1./xdata[i] for i in range(len(xdata))])
'''

#plt.figure()
#plt.Axes.fill_between(xdata,y1dat,y2dat)
regionfig,(ax) = plt.subplots()
ax.fill_between(xdatatot,y1dattot,y2dattot,where=y1dattot>y2dattot,interpolate=True,facecolor='greenyellow')
ax.fill_between(xdatapos2,2,y1dattot2,interpolate=True,facecolor='lightblue')
ax.plot([1], [1], 'ro')
#ax.plot([.1], [.1], 'o', color='orange')
ax.plot([.3], [.3], 'go')#around K 3/4
ax.plot([.5], [.5], 'bo')#    at K 2/3
ax.plot([.8], [.8], 'mo')#around K 6/11
ax.plot([1.2], [1.2], 'o', color='orange')#around K 5/11
#plt.rc('sans-serif', sans-serif=['Verdana'])
fntsz=12
#ax.text(0.2,0.1, 'I', fontsize=fntsz, ha='center')
#ax.text(-0.15,0.1, 'II', fontsize=fntsz, ha='center')
#ax.text(-0.15,-0.25, 'III', fontsize=fntsz, ha='center')
#ax.text(0.2,-0.25, 'IV', fontsize=fntsz, ha='center')
ax.text(0.3, 0.7, '1', fontsize=fntsz, ha='center')
ax.text(-.8, 0.5, '2', fontsize=fntsz, ha='center')
ax.text(-.8, -.8, '3', fontsize=fntsz, ha='center')
ax.text(0.5, -.8, '2', fontsize=fntsz, ha='center')
ax.text(1.5, 0.5, '5', fontsize=fntsz, ha='center')
ax.text(1.5, 1.5, '4', fontsize=fntsz, ha='center')
ax.text(0.5, 1.5, '5', fontsize=fntsz, ha='center')
ax.text(-.8, 1.5, '6', fontsize=fntsz, ha='center')
ax.text(-1.5, -1.5, '7', fontsize=fntsz, ha='center')
ax.text(1.5, -.8, '6', fontsize=fntsz, ha='center')
#ax.grid(True, which='both') #,linestyle='-'
ax.axvline(x=0, color='k',lw=0.5)
ax.axhline(y=0, color='k',lw=0.5)
#plt.ytick(range(-2,2))
plt.ylabel(r'$a_{12}$',fontsize=fntsz+4)
plt.xlabel(r'$a_{21}$',fontsize=fntsz+4)
plt.xlim(-2,2)
plt.ylim(-2,2)
ax.set_aspect(aspect='equal')
ax.plot([0], [0], 'o', color='grey')
#plt.legend()
plt.show()
#regionfig.savefig("a-a-graph7.pdf",bbox_inches='tight')

"""
potfig=plt.figure(2)
tempK=100;
tempC=1.6*tempK;
potx = np.arange(0,tempC,0.5)
poty = -potx*potx/2 + potx**3/3/tempK
plt.plot(potx,poty,'-')#,label='$a$ = %.1f'%(i/float(numberofas-1))
plt.ylabel('logistic potential')
plt.xlabel('population size, $n$')
plt.xlim(0,tempC)
#plt.ylim(0,4)
#plt.legend()
plt.show()
potfig.savefig("logistic-potential.pdf",bbox_inches='tight')
"""

'old'
"""
fntsz=9
ax.annotate('Moran limit', xy=(1,1), xytext=(-0.95,1.5),arrowprops=dict(facecolor='black',shrink=0.02),fontsize=fntsz)
#ax.text(-1, -0.5, 'mutualism', fontsize=fntsz)
ax.text(-1.5, -0.4, 'mutualism/symbiosis', fontsize=fntsz, ha='center')
ax.text(0.5, -0.6, 'weak\nparasitism', fontsize=fntsz, ha='center')
ax.text(-1.5, 0.5, 'weak\nparasitism', fontsize=fntsz, ha='center')
ax.text(1.5, -0.6, 'strong\nparasitism', fontsize=fntsz, ha='center')
ax.text(-1.5, 1.2, 'strong\nparasitism', fontsize=fntsz, ha='center')
ax.text(0.5, 0.3, 'weak\ncompetition\n(coexistence)', fontsize=fntsz, ha='center')
ax.text(1.5, 1.3, 'strong\ncompetition\n(competitive\nexclusion)', fontsize=fntsz, ha='center')
"""
'also old'
"""
#plt.figure()
#plt.Axes.fill_between(xdata,y1dat,y2dat)
regionfig,(ax) = plt.subplots()
ax.fill_between(xdatatot,y1dattot,y2dattot,where=y1dattot>y2dattot,interpolate=True,facecolor='greenyellow')
ax.fill_between(xdatapos2,2,y1dattot2,interpolate=True,facecolor='hotpink')
ax.plot([1], [1], 'ro')
#ax.plot([.1], [.1], 'o', color='orange')
ax.plot([.3], [.3], 'go')#around K 3/4
ax.plot([.5], [.5], 'bo')#    at K 2/3
ax.plot([.8], [.8], 'mo')#around K 6/11
ax.plot([1.2], [1.2], 'o', color='orange')#around K 5/11
#plt.rc('sans-serif', sans-serif=['Verdana'])
fntsz=10
ax.text(0.1,0.1, 'I', fontsize=fntsz)
ax.text(-0.2,0.1, 'II', fontsize=fntsz)
ax.text(-0.2,-0.2, 'III', fontsize=fntsz)
ax.text(0.1,-0.2, 'IV', fontsize=fntsz)
ax.annotate('Moran limit', xy=(1,1), xytext=(-1.7,1.7),arrowprops=dict(facecolor='black',shrink=0.00,width=0.7,headwidth=8),fontsize=fntsz+4)
ax.text(-1.0, -0.5, 'weak mutualism\n(coexistence)', fontsize=fntsz, ha='center')
#ax.text(-1.0, -0.4, 'weak mutualism\n(coexistence)', fontsize=fntsz, ha='center', fontproperties=font0)
ax.text(-1.3, -1.85, 'strong mutualism\n(pop. explosion)', fontsize=fntsz, ha='center')
ax.text(0.5, -0.8, 'weak\nparasitism\n(coexistence)', fontsize=fntsz, ha='center')
ax.text(-1.0, 0.5, 'weak parasitism\n(coexistence)', fontsize=fntsz, ha='center')
ax.text(1.5, -1.5, 'strong\nparasitism\n(domination)', fontsize=fntsz, ha='center')
ax.text(-1.0, 1.1, 'strong parasitism\n(domination)', fontsize=fntsz, ha='center')
#ax.text(0.5, 0.3, 'weak\ncompetition\n(coexistence)', fontsize=fntsz, ha='center')
ax.annotate('weak\ncompetition\n(coexistence)', xy=(.7,.2), xytext=(1.5,.4),arrowprops=dict(facecolor='black',shrink=0.00,width=.5,headwidth=5,headlength=8),fontsize=fntsz, ha='center')
ax.text(1.0, 1.5, 'strong competition\n(competitive exclusion)', fontsize=fntsz, ha='center')
#ax.grid(True, which='both') #,linestyle='-'
ax.axvline(x=0, color='k',lw=0.5)
ax.axhline(y=0, color='k',lw=0.5)
#plt.ytick(range(-2,2))
plt.ylabel(r'$a_{12}$',fontsize=fntsz+4)
plt.xlabel(r'$a_{21}$',fontsize=fntsz+4)
plt.xlim(-2,2)
plt.ylim(-2,2)
ax.set_aspect(aspect='equal')
ax.plot([0], [0], 'o', color='grey')
#plt.legend()
plt.show()
regionfig.savefig("a-a-graph6.pdf",bbox_inches='tight')
"""
"""
xdatapos2 = np.arange(1,2,0.01)
y1dattot2 = np.array([1. for i in range(len(xdatapos2))])
#y2datneg = 1./xdataneg
#y2datpos = np.array([-2 for i in range(len(xdatapos))])
#y2dattot = np.concatenate((y2datneg,y2datpos))
regionfig,(ax) = plt.subplots()
ax.fill_between(xdatatot,y1dattot,y2dattot,where=y1dattot>y2dattot,interpolate=True)
ax.fill_between(xdatapos2,2,y1dattot2,interpolate=True)
plt.ylabel(r'$a_{12}$',fontsize=14)
plt.xlabel(r'$a_{21}$',fontsize=14)
ax.axvline(x=0, color='k',lw=0.5)
ax.axhline(y=0, color='k',lw=0.5)
plt.xlim(-2,2)
plt.ylim(-2,2)
ax.set_aspect(aspect='equal')
plt.show()
#regionfig.savefig("fp-physical.pdf",bbox_inches='tight')

xdatapos3 = np.arange(0.,1,0.01)
y1dattot3 = 1./xdatapos3
xdatapos4 = np.arange(1,2,0.01)
y1dattot4 = 1./xdatapos4
regionfig,(ax) = plt.subplots()
ax.fill_between(xdatatot,y1dattot,y2dattot,where=y1dattot>y2dattot,interpolate=True, facecolor='yellowgreen')
ax.fill_between(xdatapos3,2,y1dattot3,interpolate=True, facecolor='greenyellow')
ax.fill_between(xdatapos4,y1dattot4,1,interpolate=True, facecolor='hotpink')
plt.ylabel(r'$a_{12}$',fontsize=14)
plt.xlabel(r'$a_{21}$',fontsize=14)
ax.axvline(x=0, color='k',lw=0.5)
ax.axhline(y=0, color='k',lw=0.5)
plt.xlim(-2,2)
plt.ylim(-2,2)
ax.set_aspect(aspect='equal')
plt.show()
#regionfig.savefig("fp-stable.pdf",bbox_inches='tight')
"""
