# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 18:45:31 2023

@author: moure
"""

import numpy as np
from scipy.stats import norm

def EuCallPayOff(S, K):
    if S>=K:
        return S-K
    else:
        return 0

def EuPutPayOff(S, K):
    if S<=K:
        return K-S
    else:
        return 0

def BinCallPayOff(S, K):
    if S>=K:
        return 1
    else:
        return 0

def BinPutPayOff(S, K):
    if S<=K:
        return 1
    else:
        return 0
    
def AnalyticOptionPricing(r, sigma, S0, T, payOff, K):
    d1 = (np.log(S0/K)+(r+sigma*sigma/2)*T)/(sigma*np.sqrt(T))
    d2 = d1 - sigma*np.sqrt(T)
    if payOff=="European Call":
        optValue = S0*norm.cdf(d1) - K*np.exp(-r*T)*norm.cdf(d2)
    elif payOff=="European Put":
        optValue = -S0*norm.cdf(-d1) + K*np.exp(-r*T)*norm.cdf(-d2)
    elif payOff=="Binary Call":
        optValue = np.exp(-r*T)*norm.cdf(d2)
    elif payOff=="Binary Put":
        optValue = np.exp(-r*T)*norm.cdf(-d2)
    else:
        pass
    return optValue
    
def MCOptionPricing(nbSimul, r, sigma, S0, T, payOff, timeTenor, K):
    deltat = T/timeTenor
    asset = np.zeros((nbSimul, timeTenor))
    asset[:,0] = S0
    optValue = 0
    for k in range(0,nbSimul):
        for j in range(1, timeTenor):
            asset[k][j] = asset[k][j-1]*(r*deltat+1+sigma*np.sqrt(deltat)*np.random.normal(0,1))
    if payOff=="European Call":
        temp = np.zeros(nbSimul)
        for i in range(0,nbSimul):
            temp[i] = EuCallPayOff(asset[i,timeTenor-1],K)
        optValue = np.exp(-r*T)*np.mean(temp)
    elif payOff=="European Put":
        temp = np.zeros(nbSimul)
        for i in range(0,nbSimul):
            temp[i] = EuPutPayOff(asset[i,timeTenor-1],K)
        optValue = np.exp(-r*T)*np.mean(temp)
    elif payOff=="Binary Call":
        temp = np.zeros(nbSimul)
        for i in range(0,nbSimul):
            temp[i] = BinCallPayOff(asset[i,timeTenor-1],K)
        optValue = np.exp(-r*T)*np.mean(temp)
    elif payOff=="Binary Put":
        temp = np.zeros(nbSimul)
        for i in range(0,nbSimul):
            temp[i] = BinPutPayOff(asset[i,timeTenor-1],K)
        optValue = np.exp(-r*T)*np.mean(temp)
    else:
        pass
    return optValue

def PDECallOptionPricing(assetRange, r,divYield, sigma, S0, T,timeTenor, K):
    deltaS = 2*K/assetRange
    deltaT = T/timeTenor
    timeSize = np.floor(T/deltaT).astype(int)
    assetTab = np.zeros(assetRange)
    callOptTab = np.zeros((timeSize, assetRange))
    deltaTab = np.zeros((timeSize, assetRange))
    gammaTab = np.zeros((timeSize, assetRange))
    thetaTab = np.zeros((timeSize, assetRange))
    
    #Initial fills and boundaries
    for k in range(0, assetRange):
        assetTab[k] = k*deltaS
        callOptTab[timeSize-1][k] = EuCallPayOff(assetTab[k],K)
        deltaTab[timeSize-1][k] = BinCallPayOff(assetTab[k],K)

    #Full algorithm
    for p in reversed(range(0,timeSize-1)):
        for k in reversed(range(0, assetRange)):
            if (k==0):
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k])/(deltaS)
            elif (k==(assetRange-1)):
                deltaTab[p][k] = (callOptTab[p+1][k]- callOptTab[p+1][k-1])/(deltaS)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
            else:
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k-1])/(2*(deltaS))
                gammaTab[p][k] = (callOptTab[p+1][k+1]+ callOptTab[p+1][k-1] - 2*callOptTab[p+1][k])/(deltaS**2)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
    proxy = int(np.ceil(S0/2))
    return callOptTab[0,:][proxy]

def PDEPutOptionPricing(assetRange, r,divYield, sigma, S0, T,timeTenor, K):
    deltaS = 2*K/assetRange
    deltaT = T/timeTenor
    timeSize = np.floor(T/deltaT).astype(int)
    assetTab = np.zeros(assetRange)
    callOptTab = np.zeros((timeSize, assetRange))
    deltaTab = np.zeros((timeSize, assetRange))
    gammaTab = np.zeros((timeSize, assetRange))
    thetaTab = np.zeros((timeSize, assetRange))
    
    #Initial fills and boundaries
    for k in range(0, assetRange):
        assetTab[k] = k*deltaS
        callOptTab[timeSize-1][k] = EuPutPayOff(assetTab[k],K)
        deltaTab[timeSize-1][k] = BinPutPayOff(assetTab[k],K)

    #Full algorithm
    for p in reversed(range(0,timeSize-1)):
        for k in reversed(range(0, assetRange)):
            if (k==0):
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k])/(deltaS)
            elif (k==(assetRange-1)):
                deltaTab[p][k] = (callOptTab[p+1][k]- callOptTab[p+1][k-1])/(deltaS)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
            else:
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k-1])/(2*(deltaS))
                gammaTab[p][k] = (callOptTab[p+1][k+1]+ callOptTab[p+1][k-1] - 2*callOptTab[p+1][k])/(deltaS**2)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
    proxy = int(np.ceil(S0/2))
    return callOptTab[0,:][proxy]

def PDEBinCallOptionPricing(assetRange, r,divYield, sigma, S0, T,timeTenor, K):
    deltaS = 2*K/assetRange
    deltaT = T/timeTenor
    timeSize = np.floor(T/deltaT).astype(int)
    assetTab = np.zeros(assetRange)
    callOptTab = np.zeros((timeSize, assetRange))
    deltaTab = np.zeros((timeSize, assetRange))
    gammaTab = np.zeros((timeSize, assetRange))
    thetaTab = np.zeros((timeSize, assetRange))
    
    #Initial fills and boundaries
    for k in range(0, assetRange):
        assetTab[k] = k*deltaS
        callOptTab[timeSize-1][k] = BinCallPayOff(assetTab[k],K)

    #Full algorithm
    for p in reversed(range(0,timeSize-1)):
        for k in reversed(range(0, assetRange)):
            if (k==0):
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k])/(deltaS)
            elif (k==(assetRange-1)):
                deltaTab[p][k] = (callOptTab[p+1][k]- callOptTab[p+1][k-1])/(deltaS)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
            else:
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k-1])/(2*(deltaS))
                gammaTab[p][k] = (callOptTab[p+1][k+1]+ callOptTab[p+1][k-1] - 2*callOptTab[p+1][k])/(deltaS**2)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
    proxy = int(np.ceil(S0/2))
    return callOptTab[0,:][proxy]

def PDEBinPutOptionPricing(assetRange, r,divYield, sigma, S0, T,timeTenor, K):
    deltaS = 2*K/assetRange
    deltaT = T/timeTenor
    timeSize = np.floor(T/deltaT).astype(int)
    assetTab = np.zeros(assetRange)
    callOptTab = np.zeros((timeSize, assetRange))
    deltaTab = np.zeros((timeSize, assetRange))
    gammaTab = np.zeros((timeSize, assetRange))
    thetaTab = np.zeros((timeSize, assetRange))
    
    #Initial fills and boundaries
    for k in range(0, assetRange):
        assetTab[k] = k*deltaS
        callOptTab[timeSize-1][k] = BinPutPayOff(assetTab[k],K)

    #Full algorithm
    for p in reversed(range(0,timeSize-1)):
        for k in reversed(range(0, assetRange)):
            if (k==0):
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k])/(deltaS)
            elif (k==(assetRange-1)):
                deltaTab[p][k] = (callOptTab[p+1][k]- callOptTab[p+1][k-1])/(deltaS)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
            else:
                deltaTab[p][k] = (callOptTab[p+1][k+1]- callOptTab[p+1][k-1])/(2*(deltaS))
                gammaTab[p][k] = (callOptTab[p+1][k+1]+ callOptTab[p+1][k-1] - 2*callOptTab[p+1][k])/(deltaS**2)
                thetaTab[p][k] = r*callOptTab[p+1][k] - (r-divYield)*assetTab[k]*deltaTab[p][k]-sigma**2*assetTab[k]**2*gammaTab[p][k]/2
                callOptTab[p][k] = callOptTab[p+1][k] - thetaTab[p][k]*deltaT
    proxy = int(np.ceil(S0/2))
    return callOptTab[0,:][proxy]


print(MCOptionPricing(1000, 0.02, 0.1, 100,1,'European Call',10,100))
print(MCOptionPricing(1000, 0.02, 0.1, 100,1,'European Put',10,100))
print(MCOptionPricing(1000, 0.02, 0.1, 100,1,'Binary Call',10,100))
print(MCOptionPricing(1000, 0.02, 0.1, 100,1,'Binary Put',10,100))
print(AnalyticOptionPricing(0.02, 0.1, 100,1,'European Call',100))
print(AnalyticOptionPricing(0.02, 0.1, 100,1,'European Put',100))
print(AnalyticOptionPricing(0.02, 0.1, 100,1,'Binary Call',100))
print(AnalyticOptionPricing(0.02, 0.1, 100,1,'Binary Put',100))
print(PDECallOptionPricing(100, 0.02,0, 0.1, 100,1,100,100))
print(PDEPutOptionPricing(100, 0.02,0, 0.1, 100,1,100,100))
print(PDEBinCallOptionPricing(100, 0.02,0, 0.1, 100,1,100,100))
print(PDEBinPutOptionPricing(100, 0.02,0, 0.1, 100,1,100,100))
