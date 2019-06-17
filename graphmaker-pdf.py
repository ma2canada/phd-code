import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as spint
from scipy.special import poch

def birthrate(n,d,q,K):
        return n*(1+d/2)-q*n*n/K
def deathrate(n,d,q,K):
        return n*d/2+(1-q)*n*n/K

def pdfFP_gaussian(x, stoch, cap, delta=1.):
    return np.exp((-(cap-x)**2)/(cap*(2+delta-2*stoch)))
def pdfFP_gaussian_normalized(x, stoch, cap, delta=1.):
    return np.exp((-(cap-x)**2)/(cap*(2+delta-2*stoch)))/np.sqrt(np.pi*cap*(2+delta-2*stoch))
def mteFP_gaussian(cap,stoch,alsocap,bigN,delta=1.):
    return 1./(deathrate(1,delta,stoch,cap)*pdfFP_gaussian_normalized(1,stoch,cap,delta))

def pdfFP_full(x, stoch, cap, delta=1.):
#    temp = np.exp(-2*cap/(1-2*stoch))*(1/cap)*((cap*(1-2*stoch)+cap*(1+delta))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))
#    return np.exp(-2*x/(1-2*stoch))*(1/x)*((x*(1-2*stoch)+cap*(1+delta))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))/temp
    return np.exp(-2*x/(1-2*stoch)+2*cap/(1-2*stoch))*(cap/x)*(((x*(1-2*stoch)+cap*(1+delta))/((cap*(1-2*stoch)+cap*(1+delta))))**(1+(2*cap*(2+delta-2*stoch))/(1-2*stoch)**2))
def normalization_constant(stoch, cap, delta=1.):
    return spint.quad(lambda x: pdfFP_full(x,stoch,cap,delta), 0, 2.*cap)[0]
def pdfFP_full_normalized(x, stoch, cap, delta=1.):
    return pdfFP_full(x, stoch, cap, delta)/normalization_constant(stoch,cap,delta)

def grand(y, x, stoch, cap, delta):
    return (2*x - 2*x*x/cap - (1+delta) - 2*(1-2*stoch)*x/cap)/(x*(1+delta)+(1-2*stoch)*x*x/cap)
def pdfFP_full_numericarray_unnormalized(stoch, cap, delta):
    return spint.odeint(grand,[0.00001],range(1,2*cap+1),args=(stoch,cap,delta))
def pdfFP_full_numericarray(stoch, cap, delta=1.):
    return pdfFP_full_numericarray_unnormalized(stoch, cap, delta)/max(pdfFP_full_numericarray_unnormalized(stoch, cap, delta))

def pdf_smalln(x, stoch, cap, delta=1.):
    lam = 2./(2.+delta) # < 1.
#    lam = delta/(2.+delta) # < 1. #I think this is the wrong one
#    eff = lam**(x) #this is not the theory, but gives the right direction
    eff = lam**(-x)
    return eff/x
def pdf_smalln_recursive_unnormalized(x, stoch, cap, delta=1.):
    P0 = 1.E-15 #an arbitrary P0 that can probably be tweaked to match boundaries
    if x==1:
        return P0
    else:
#        return x*(1.+delta/2.-x*stoch/cap)/(x+1)/(delta/2+(1-stoch)*(x+1)/cap)*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)
        return (birthrate(x-1,delta,stoch,cap))*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)/deathrate(x,delta,stoch,cap)
#        return (birthrate(x-1,delta,stoch,cap)+deathrate(x-1,delta,stoch,cap))*pdf_smalln_recursive_unnormalized(x-1, stoch, cap, delta)/deathrate(x,delta,stoch,cap)
def pdf_smalln_recursive_normalizer(stoch, cap, delta=1.):
    '''
    if(type(cap)==float or type(cap)==int):
        normalizer = 0.
        for ex in range(1,int(cap*(1+delta)/stoch)):
            normalizer+=pdf_smalln_recursive_unnormalized(ex, stoch, cap, delta)
#    elif(type(cap)==list or type(cap)==np.ndarray):
    elif(len(cap)>1):
        normalizer = [0. for j in cap];
        for k,pacity in enumerate(cap):
            for ex in range(1,int(pacity*(1+delta)/stoch)):
                normalizer[k]+=pdf_smalln_recursive_unnormalized(ex, stoch, pacity, delta)
    else:
        print("I can't let you do that, Dave")
    return normalizer
    '''
    normalizer = 0.
    for ex in range(1,int(cap*(1+delta)/stoch)):
        normalizer+=pdf_smalln_recursive_unnormalized(ex, stoch, cap, delta)
    return normalizer
def pdf_smalln_recursive(x, stoch, cap, delta=1.):
    return pdf_smalln_recursive_unnormalized(x, stoch, cap, delta)/pdf_smalln_recursive_normalizer(stoch,cap,delta)
def pdf_smalln_recursive_list(x, stoch, cap, delta=1.):
    temp=[0. for i in x];
    for j,n in enumerate(x):
        temp[j]=pdf_smalln_recursive_unnormalized(n, stoch, cap, delta)/pdf_smalln_recursive_normalizer(stoch,cap,delta)
    return temp
def mte_smalln_recursive(stoch,cap,delta=1.):
        return pdf_smalln_recursive_normalizer(stoch,cap,delta)/deathrate(1,delta,stoch,cap)
def mte_smalln_recursive_list(stoch,cap,delta=1.):
        temp=[0. for i in cap];
        for j,pacity in enumerate(cap):
            temp[j]=1./(deathrate(1,delta,stoch,pacity)*pdf_smalln_recursive(1,stoch,pacity,delta))
#            temp[j]=pdf_smalln_recursive_normalizer(stoch,pacity,delta)/deathrate(1,delta,stoch,pacity)
        return temp

def pdf_smalln_better_unnormalized(x, stoch, cap, delta=1.):#this creates problems as it divides ~0 by ~0 or something like that
    P1 = 1.;
    return P1*(1-stoch)/(1-2*stoch)/x*poch(1+cap*(1+delta)/(1-2*stoch),x-1)/poch(1+(1-stoch+cap*delta/2.)/(1-stoch),x-1)*((1-2*stoch)/(1-stoch))**x
def pdf_smalln_better_P1(stoch, cap, delta=1.):
    normalizer = 0.
    for ex in range(1,int(cap*(1+delta)/stoch)):
        normalizer+=pdf_smalln_better_unnormalized(ex, stoch, cap, delta)
    return 1./normalizer
def pdf_smalln_better(x, stoch, cap, delta=1.):
    return pdf_smalln_better_unnormalized(x, stoch, cap, delta)*pdf_smalln_better_P1(stoch,cap,delta)
def pdf_smalln_better_recursive_unnormalized(x, stoch, cap, delta=1.):
    P0 = 1.E-15 #an arbitrary P0 that can probably be tweaked to match boundaries
    if x==1:
        return P0
    else:
        return (birthrate(x-1,delta,stoch,cap)+deathrate(x-1,delta,stoch,cap))/deathrate(x,delta,stoch,cap)*pdf_smalln_better_recursive_unnormalized(x-1, stoch, cap, delta)
def pdf_smalln_better_recursive_normalizer(stoch, cap, delta=1.):
    normalizer = 0.
    for ex in range(1,int(cap*(1+delta)/stoch)):
        normalizer+=pdf_smalln_better_recursive_unnormalized(ex, stoch, cap, delta)
    return normalizer
def pdf_smalln_better_recursive(x, stoch, cap, delta=1.):
    return pdf_smalln_better_recursive_unnormalized(x, stoch, cap, delta)/pdf_smalln_better_recursive_normalizer(stoch,cap,delta)
def mte_smalln_better(stoch,cap,delta=1.):
    return 1./(deathrate(1,delta,stoch,cap)*pdf_smalln_better_P1(stoch,cap,delta))
def mte_smalln_better_recursive(stoch,cap,delta=1.):
    return 1./(deathrate(1,delta,stoch,cap)*pdf_smalln_better_recursive_normalizer(stoch,cap,delta))
#    return 1./(deathrate(1,delta,stoch,cap)*pdf_smalln_better_recursive(1,stoch,cap,delta))
def mte_smalln_better_recursive_list(stoch,cap,delta=1.):
    temp = [1. for i in cap];
    for j,pacity in enumerate(cap):
        temp[j] = 1./(deathrate(1,delta,stoch,pacity)*pdf_smalln_better_recursive(1,stoch,pacity,delta))
#        temp[j] = pdf_smalln_better_recursive_normalizer(stoch,pacity,delta)*1./(deathrate(1,delta,stoch,pacity))
    return temp

#------------WKB realspace------------------
def WKB_RS(n, stoch, cap, N, delta):
    action =  (cap/N/n/(1-stoch))*deathrate(n,delta,stoch,cap)*np.log(cap*deathrate(n,delta,stoch,cap)/n) + (cap/N/stoch/n)*birthrate(n,delta,stoch,cap)*np.log(cap*birthrate(n,delta,stoch,cap)/n) - (cap/N/n)*(deathrate(n,delta,stoch,cap)/(1-stoch) + birthrate(n,delta,stoch,cap)/stoch)
    action1 = (1/2)*np.log( deathrate(n, delta, stoch, cap)*birthrate(n, delta, stoch, cap) )
    secDerivAction = N*( (delta/2+(1-stoch)*2*n/cap)/deathrate(n, delta, stoch, cap) - (1+delta/2-stoch*2*n/cap)/birthrate(n, delta, stoch, cap) ) #ie. N*(d'/d - b'/b)
    return action, action1, secDerivAction
def distributionWKB_RS(n, fixPoints, stoch, cap, N, delta):
    act, act1, secDerivAct = WKB_RS(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = WKB_RS(fixPoints, stoch, cap, N, delta)
    #constant term in the probability distribution found from WKB
    distribution = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(-N*(act-actF) - (act1-act1F))
    return distribution
def mteWKB(fixedPoints, stoch, cap, N, delta): #MAB
    return 1/(deathrate(1, delta, stoch, cap)*distributionWKB_RS(1, fixedPoints, stoch, cap, N, delta))

#------------_ACTUALLY_ FP WKB------------------ #MAB
def actionFP(n, stoch, cap, N, delta): #MAB
    action = 2.*(1/N)*(-1/(1-2*stoch)**2)*( cap*(2+delta-2*stoch)*np.log( (cap*(1+delta)+(1-2*stoch)*n)/(cap*(2+delta-2*stoch)) ) - (1-2*stoch)*(n-cap))
    action1 = 0.0 #this is not necessarily true, but likely doesn't matter
    secDerivAction = (2*N/cap)/(2+delta-2*stoch)
    return action, action1, secDerivAction
def statDistributionFP_quasi(n, fixPoints, stoch, cap, N, delta): #MAB
    act, act1, secDerivAct = actionFP(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = actionFP(fixPoints, stoch, cap, N, delta)
    constant = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(N*actF + act1F) #constant term in the probability distribution found from WKB
    distribution = constant*np.exp(-N*act - act1)
    return distribution
def mteFP_quasi(fixedPoints, stoch, cap, N, delta): #MAB
    return 1/(deathrate(1, delta, stoch, cap)*statDistributionFP_quasi(1, fixedPoints, stoch, cap, N, delta))

def rates(stoch, cap, delta, maxSum):
    """
    Function that calculates out outputs different ratios of the product of death rates and product of birth rates.

    Args:
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity
        maxSum(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        Q: Some weird ratio of the rates, product(birth(i-1)...birth(1)/death(i)...death(1)).
        S: Some weird ratio of the rates, product(death(i)...death(1)/birth(i)...birth(1)).
    """
    Q = np.array([deathrate(1, delta, stoch, cap)**(-1)])
    S = np.array([deathrate(1, delta, stoch, cap)*birthrate(1, delta, stoch, cap)**(-1)])
    for i in range(2,maxSum):
        Q = np.append(Q,Q[-1]*birthrate(i-1, delta, stoch, cap)*deathrate(i, delta, stoch, cap)**(-1))
        if i <= cap:
            S = np.append(S,S[-1]*deathrate(i, delta, stoch, cap)*birthrate(i, delta, stoch, cap)**(-1))
    """
    Q1 = [deathrate(1, delta, stoch, cap)**(-1)]
    S1 = [deathrate(1, delta, stoch, cap)*birthrate(1, delta, stoch, cap)**(-1)]

    for i in range(2,maxSum):
        Q1.append(Q1[-1]*birthrate(i-1, delta, stoch, cap)*deathrate(i, delta, stoch, cap)**(-1))
        if i <= cap:
            S1.append(S1[-1]*deathrate(i, delta, stoch, cap)*birthrate(i, delta, stoch, cap)**(-1))
    """
    return Q, S
def mteSum1D(fixedPoints, stoch, cap, maxSum, delta):
    """
    Function that calculates the Mean Time Extinction from different rate equations.

    Args:
        stoch(array or int): The stochastic variable in our equations
        cap(array or int): The capacity
        maxSum(array): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        mte: Mean Time Extinction.
    """
    mte = []
    if np.size(stoch)>1:
        for i, stoch_i in enumerate(stoch):
            q, s = rates(stoch_i, cap, delta, np.int(maxSum[i]))
            mte.append(np.sum(q))
            for j, ratio in enumerate(s):
                mte[i] += ratio*np.sum(q[j+1:-1])
    elif np.size(cap)>1:
        for i, cap_i in enumerate(cap):
            q, s = rates(stoch, cap_i, delta, np.int(maxSum[i]))
            mte.append(np.sum(q))
            for j, ratio in enumerate(s):
                mte[i] += ratio*np.sum(q[j+1:-1])
    elif np.size(delta)>1:
        for i, delta_i  in enumerate(delta):
            q, s = rates(stoch, cap, delta_i, np.int(maxSum[i]))
            mte.append(np.sum(q))
            for j, ratio in enumerate(s):
                mte[i] += ratio*np.sum(q[j+1:-1])
    return mte
'''
def actionWKB(stoch, cap, N, delta):
    action = (1/N)*( n*np.log( deathrate(n, delta, stoch, cap)/birthrate(n, delta, stoch, cap)) + (cap*delta/(2*(1-stoch))) * np.log((deathrate(n, delta, stoch, cap))/n) + (cap*(1+delta/2)/stoch) * np.log(birthrate(n, delta, stoch, cap))/n)
    action1 = (1/2)*np.log( deathrate(n, delta, stoch, cap)*birthrate(n, delta, stoch, cap) )
    secDerivAction = n*N/cap*( (1-stoch)/deathrate(n, delta, stoch, cap) + stoch/birthrate(n, delta, stoch, cap) )
    return action, action1, secDerivAction
def statDistributionWKB(n, fixPoints, stoch, cap, N, delta):
    act, act1, secDerivAct = actionWKB(n, stoch, cap, N, delta)
    actF, act1F, secDerivActF = actionWKB(fixPoints, stoch, cap, N, delta)
    constant = np.sqrt(secDerivActF/(2*np.pi*N))*np.exp(N*actF + act1F) #constant term in the probability distribution found from WKB
    distribution = constant*np.exp(-N*act - act1)
    return distribution
'''

'#############################################################################'
#------------Algorithm, as per Jeremy----------------

def statDistributionAlgo(fixedPoints, stoch, cap, N, delta, std=10):
    """
    Algorithm that calculates the quasistationary conditional probability distribution function

    Args:
        fixPoints(int): The fixed points of our equations (mean of intial guess for distribution).
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity of our population.
        N(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The quasistationary conditional probability distribution function
    """
    def birthrate(n, delta, stoch, cap):
        return (1+delta/2)*n-stoch*n**2/cap
    def deathrate(n, delta, stoch, cap):
        return delta*n/2+(1-stoch)*n**2/cap
    def gaussian(x, mean, std):
        return np.sqrt(2*np.pi*std**2)**(-1) * np.exp( -( (x-mean)/(np.sqrt(2)*std) )**2 )
    popDist = np.asarray([gaussian(n, cap, std) for n in range(0,N+1)])
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])

    #popDist = np.insert(popDist,0,0)
    #birthArray = np.insert(birthArray,0,0)
    #deathArray = np.insert(deathArray,0,0)
    popDist = np.append(popDist,0)
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)

    for i in range(0,500):
        statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        popDist[1:-1] = np.abs(statDist)#/sum(statDist))

    print("is it normalized? ",sum(popDist[1:-1]))
    return popDist[1:-1]

#------------my tortured version of the Algorithm----------------
def statDistributionAlgo2(fixedPoints, stoch, cap, N, delta, std=10):
    """
    Algorithm that calculates the quasistationary conditional probability distribution function

    Args:
        fixedPoints(int): The fixed points of our equations (mean of intial guess for distribution).
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity of our population.
        N(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The quasistationary conditional probability distribution function
    """
    def birthrate(n,d,q,K):
        return n*(1+d/2)-q*n*n/K
    def deathrate(n,d,q,K):
        return n*d/2+(1-q)*n*n/K
    #popDist = np.asarray([gaussian(n, cap, std) for n in range(0,N+1)])
    popDist = np.asarray([1./N for n in range(0,N+1)]) #from n=0 to n=N
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])
    popDist = np.insert(popDist,0,0) #from n=-1 to n=N
    birthArray = np.insert(birthArray,0,0)
    deathArray = np.insert(deathArray,0,0)
    popDist = np.append(popDist,0) #from n=-1 to n=N+1
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)
    """
    distDiff = 100
    distDiffPrev = 0
    while distDiff != distDiffPrev:
        print distDiff
        distDiffPrev = distDiff
        statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        distDiff = np.sum(np.abs(statDist/max(statDist)-popDist[1:-1]))
        popDist[1:-1] = np.abs(statDist/max(statDist))
    """
    #print(deathArray)
#    statDist = np.asarray([np.exp(-N) for i in range(1,N+1)]) #from n=1 to n=N #MAB
    statDist = np.asarray([1/N for i in range(1,N+1)]) #from n=1 to n=N #MAB
#    print("before normalizing: ",sum(statDist))
#    print("P^C(1) = ",statDist[0])
    for i in range(0,500):
        ###statDist goes from n=0 to n=N as written so far, for Jeremy
        #statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        ###statDist goes from n=1 to n=N as written so far, for #MAB
        #statDist[1:] = (-birthArray[1+2-2:-1-2]*popDist[1:-3] + (birthArray[2:-2]+deathArray[2:-2])*popDist[2:-2] - deathArray[2]*popDist[2:-2]*popDist[2]) / (deathArray[3:-1])#MAB
        for j in range(1,len(statDist)):#MAB
#            statDist[j] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
            popDist[j+2] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
        #popDist[1:-1] = np.abs(statDist/max(statDist))
#        if(i%100==0): print("before normalizing: ",sum(statDist))
#        popDist[2:-1] = np.abs(statDist/sum(statDist))#MAB
#        statDist[0]=popDist[0+2];#MAB
#        statDist = popDist[2:-1]#MAB
        popDist = popDist/sum(popDist)#MAB
#        if(i%100==0): print("P^C(1) = ",statDist[0])
#    print(len(statDist)," vs ",len(popDist))
#    print(statDist[len(statDist)]," vs ",popDist[len(statDist)+2])
    print("is it normalized? ",sum(popDist[2:-1]))
    print("and the other details? ",sum(popDist),popDist[0:2+1],popDist[-2:])
#    return popDist[1:-1]
    return popDist[2:-1]#MAB - since we want from n=1 
#    return popDist[2:]#MAB - this gives from n=0 to n=N+1 
    '''
    #print(deathArray)
    statDist = np.asarray([1/N for i in range(1,N+1)])#MAB
    for i in range(0,5000):
        ###statDist goes from n=0 to n=N as written so far, for Jeremy
        #statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        ###statDist goes from n=1 to n=N as written so far, for #MAB
#        statDist = (-birthArray[:-3]*popDist[:-3] + (birthArray[1:-2]+deathArray[1:-2])*popDist[1:-2] - deathArray[2]*popDist[1:-2]*popDist[2]) / (deathArray[2:-1])#MAB
        for j in range(len(statDist)):#MAB
            #statDist[j] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
            popDist[j+2] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
        #print(len(statDist)," vs ",len(popDist))
        #popDist[1:-1] = np.abs(statDist/max(statDist))
        #popDist[2:-1] = np.abs(statDist/sum(statDist))#MAB
        popDist = popDist/sum(popDist)#MAB
        #print(sum(statDist))
#    return popDist[1:-1]
    return popDist[2:-1]#MAB - since we want from n=1 
#    return popDist[2:]#MAB - this gives from n=0 to n=N+1 
    '''

#------------my tortured version of the Algorithm - again----------------
def statDistributionAlgo3(fixedPoints, stoch, cap, N, delta, std=10):
    def birthrate(n,d,q,K):
        return n*(1+d/2)-q*n*n/K
    def deathrate(n,d,q,K):
        return n*d/2+(1-q)*n*n/K
    popDist = np.asarray([1./N for n in range(0,N+1)]) #from n=0 to n=N
    birthArray = np.asarray([birthrate(n, delta, stoch, cap) for n in range(0,N+1)])
    deathArray = np.asarray([deathrate(n, delta, stoch, cap) for n in range(0,N+1)])
    popDist = np.insert(popDist,0,0) #from n=-1 to n=N
    birthArray = np.insert(birthArray,0,0)
    deathArray = np.insert(deathArray,0,0)
    popDist = np.append(popDist,0) #from n=-1 to n=N+1
    birthArray = np.append(birthArray,0)
    deathArray = np.append(deathArray,0)
    statDist = np.asarray([1/N for i in range(1,N+1)]) #from n=1 to n=N #MAB
    for i in range(0,10000):
        ###statDist goes from n=0 to n=N as written so far, for Jeremy
        #statDist = (birthArray[:-2]*popDist[:-2] + deathArray[2:]*popDist[2:]) / (birthArray[1:-1] + deathArray[1:-1] - popDist[2]*deathArray[2])
        ###statDist goes from n=1 to n=N as written so far, for #MAB
        statDist[1:] = (-birthArray[1+2-2:-1-2]*popDist[1:-3] + (birthArray[2:-2]+deathArray[2:-2])*popDist[2:-2] - deathArray[2]*popDist[2:-2]*popDist[2]) / (deathArray[3:-1])#MAB
        #for j in range(1,len(statDist)):#MAB
#            statDist[j] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
            #popDist[j+2] = (-birthArray[j]*popDist[j] + (birthArray[j+1]+deathArray[j+1])*popDist[j+1] - deathArray[2]*popDist[j+1]*popDist[2]) / (deathArray[j+2])#MAB
        popDist[2:-1] = np.abs(statDist/sum(statDist))#MAB
        statDist = popDist[2:-1]#MAB
        #popDist = popDist/sum(popDist)#MAB
    print("is THIS normalized? ",sum(popDist[2:-1]))
    print("aND The other details? ",sum(popDist),popDist[0:2+1],popDist[-2:])
    return popDist[2:-1]#MAB - since we want from n=1 


def quasiStatDist(fixedPoints, stoch, cap, N, delta, std=10):
    """
    Algorithm that calculates the quasistationary conditional probability distribution function at population 1 for multiple K
    Args:
        fixedPoints(int): The fixed points of our equations (mean of intial guess for distribution).
        stoch(int): The stochastic variable in our equations
        cap(int): The capacity of our population.
        N(int): The maximal population size (determined when the birthrate is 0)
        delta: The coefficient which determines the stochasticity of the population.
    Returns:
        The quasistationary conditional probability distribution function
    """
    dist1 = np.asarray([])
    for i, capacity in enumerate(cap):
        dist = statDistributionAlgo(fixedPoints[i], stoch, capacity, N[i], delta, std)
        dist1 = np.append(dist1,dist[1])
        print(i)
    return dist1

def myQSDtime(fixedPoints, stoch, cap, N, delta, std=10):
    return 1./quasiStatDist(fixedPoints, stoch, cap, N, delta)/(delta+(1-stoch)/cap)
'#############################################################################'


makeK = 25;
makeq = 0.1;
maked = 0.01;
buff = 2;

''' some probabilities (unnormalized?)'''
"""
potfig=plt.figure(1)
'''x_array = [i for i in range(1,buff*makeK+1)];
#pdfFP_gaussian_array = [pdfFP_gaussian(i,makeq,makeK,maked) for i in range(1,101)];
#pdfFP_full_array = [pdfFP_full(i,makeq,makeK,maked) for i in range(1,101)];
pdfFP_gaussian_array = [pdfFP_gaussian_normalized(i,makeq,makeK,maked) for i in range(1,buff*makeK+1)];
pdfFP_full_array = [pdfFP_full_normalized(i,makeq,makeK,maked) for i in range(1,buff*makeK+1)];
pdfFP_WKB_array = [statDistributionFP_quasi(i,makeK,makeq,makeK,makeK,maked) for i in range(1,buff*makeK+1)];
#pdf_smalln_array = [pdf_smalln(i,makeq,makeK,maked) for i in range(1,101)];
pdf_smalln_array2 = [pdf_smalln_recursive(i,makeq,makeK,maked) for i in range(1,buff*makeK+1)];
pdf_smalln_better_array = [pdf_smalln_better(i,makeq,makeK,maked) for i in range(1,buff*makeK+1)];'''
x_array = []; pdfFP_gaussian_array = []; pdfFP_full_array = []; pdf_smalln_array = []; pdf_smalln_array2 = []; pdf_smalln_better_array = []; pdf_smalln_better_array2 = []; pdfFP_WKB_array = []; pdf_WKB_array = [];
for i in range(1,buff*makeK+1):
    x_array.append(i)
    pdfFP_gaussian_array.append(pdfFP_gaussian_normalized(i,makeq,makeK,maked))
    pdfFP_full_array.append(pdfFP_full_normalized(i,makeq,makeK,maked))
    pdf_smalln_array.append(pdf_smalln(i,makeq,makeK,maked))
    pdf_smalln_array2.append(pdf_smalln_recursive(i,makeq,makeK,maked))
    pdf_smalln_better_array.append(pdf_smalln_better(i,makeq,makeK,maked))
    pdf_smalln_better_array2.append(pdf_smalln_better_recursive(i,makeq,makeK,maked))
    pdfFP_WKB_array.append(statDistributionFP_quasi(i,makeK,makeq,makeK,makeK,maked))
    pdf_WKB_array.append(distributionWKB_RS(i,makeK,makeq,makeK,makeK,maked))

plt.semilogy(x_array,pdfFP_gaussian_array,':',label='FP Gaussian')
plt.semilogy(x_array,pdfFP_full_array,'--',label='FP numeric')
plt.semilogy(x_array,pdfFP_WKB_array,':',label='FP WKB')
plt.semilogy(x_array,pdf_WKB_array,':',label='WKB algorithm')
#plt.semilogy(x_array,pdf_smalln_array,'-',label='smalln1')#WRONG
plt.semilogy(x_array,pdf_smalln_array2,'-.',label='small $n$')
#plt.semilogy(x_array,pdf_smalln_better_array,'-',label='better smalln')#DOESN'T PLOT
#plt.semilogy(x_array,pdf_smalln_better_array2,'-',label='better smalln2')#WRONG
#plt.semilogy(x_array,statDistributionAlgo(makeK,makeq,makeK,buff*makeK,maked),':',label='QSD - Jem')
plt.semilogy(x_array,statDistributionAlgo2(makeK,makeq,makeK,buff*makeK,maked),'-',label='QSD algorithm')
#plt.semilogy(x_array,statDistributionAlgo3(makeK,makeq,makeK,buff*makeK,maked),'-',label='QSD - me3')
#plt.semilogy(x_array,pdfFP_full_numericarray(makeq,makeK,maked),'-')
plt.ylabel('stationary conditional probability distribution function')
plt.xlabel('population size $n$')
plt.xlim(0,101)
plt.ylim(1E-20,1)
plt.legend()
plt.show()
#potfig.savefig("Figure4our.pdf",bbox_inches='tight')
"""

"""
maximum = buff*makeK;
population = np.linspace(1,maximum,maximum)#all other arguments are single valued

pdfplot=plt.figure(2)
plt.semilogy(population,pdfFP_gaussian_normalized(population,makeq,makeK,maked),'-',label='gauss')
plt.semilogy(population,statDistributionFP_quasi(population,makeK,makeq,makeK,makeK,maked),'-',label='FP WKB')
plt.semilogy(population,pdf_smalln_recursive_list(population,makeq,makeK,maked),'-',label='smalln2')
plt.ylabel('pdf')
plt.xlabel('$n$')
plt.xlim(0,101)
plt.ylim(1E-20,1)
plt.legend()
plt.show()
"""

capacity = makeK;
cap = np.linspace(1.0, capacity, num=capacity)

mteplot=plt.figure(3)
plt.semilogy(mteFP_gaussian(cap,makeq,cap,cap,maked),'--',label='FP Gaussian')
#plt.semilogy(mteFP_quasi(cap,makeq,cap,cap,maked),'-',label='FP WKB')
plt.semilogy(mte_smalln_recursive_list(makeq,cap,maked),':',label='small $n$')
#plt.semilogy(mte_smalln_better(makeq,cap,maked),'-',label='better smalln')
#plt.semilogy(mte_smalln_better_recursive_list(makeq,cap,maked),'-',label='better smalln2')
plt.semilogy(mteWKB(cap,makeq,cap,cap,maked),'-.',label='WKB')
#plt.semilogy([mteSum1D(k, makeq, k, 2*k, maked) for k in cap],'-',label='1D sum')
#plt.semilogy([myQSDtime(k, makeq, k, 2*k, maked) for k in cap],'-',label='QSD (~1D sum)')
#plt.semilogy(myQSDtime(cap, makeq, cap, 2*cap, maked),'-',label='QSD (~1D sum)')
plt.semilogy(np.exp(cap)/cap,'-',label='$\delta,q=0$')
plt.ylabel('$\tau_e$')
plt.xlabel('carrying capacity $K$')
plt.xlim(0,makeK)
#plt.ylim(1E0,1E30)
plt.legend()
plt.show()

#k=1
#print(mteSum1D(k, makeq, k, 2*k, maked))
#print(myQSDtime(k, makeq, k, 2*k, maked))

