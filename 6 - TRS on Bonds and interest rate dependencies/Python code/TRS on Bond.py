# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 11:52:30 2023
@author: moureaux pierre
"""

import numpy as np
from abc import ABC, abstractmethod
from numpy import linalg as LA

"""The abstract finite difference class"""
class FiniteDifferences(object):

    def __init__(self, r0, T, sigma, alpha, beta,rmin, rmax, M, N):
        self.r0 = r0
        self.T = T
        self.sigma = sigma
        self.alpha = alpha
        self.beta = beta
        self.rmin = rmin
        self.rmax = rmax
        self.M, self.N = int(M), int(N)
        self.boundary_conds = np.linspace(rmin, rmax, self.M+1)
        self.dr = (rmax - rmin) / float(self.M)
        if self.dr > (self.sigma**2)/(LA.norm(self.alpha*(self.beta-self.boundary_conds), ord=np.inf)):
            self.dr = (self.sigma**2)/(LA.norm(self.alpha*(self.beta-self.boundary_conds), ord=np.inf))
        self.dt = T / float(self.N)
        if self.dt > (self.dr**2)/(self.sigma**2):
            self.dt = (self.dr**2)/(self.sigma**2)
        self.j_values = np.arange(self.N)
        self.grid = np.zeros(shape=(self.M+1, self.N+1))
    
    @abstractmethod
    def _setup_boundary_conditions_(self):
        pass
    
    def _setup_coefficients_(self):
        self.a = self.dt/2*((self.sigma/self.dr)**2-self.alpha*(self.beta-self.boundary_conds)/self.dr)
        self.b = -(self.sigma**2)*self.dt/(self.dr**2)-self.boundary_conds*self.dt+1
        self.c = self.dt/2*((self.sigma/self.dr)**2+self.alpha*(self.beta-self.boundary_conds)/self.dr)
        self.a0 = self.dt*(-((self.sigma/self.dr)**2)/3-1/self.dr*self.alpha*(self.beta-self.boundary_conds)-self.boundary_conds)+1
        self.b0 = self.dt/self.dr*((self.sigma**2)/(2*self.dr)+self.alpha*(self.beta-self.boundary_conds))
        self.c0 = -self.dt*(self.sigma**2)/(6*(self.dr**2))
        self.aM = self.dt*(((self.sigma/self.dr)**2)/3+1/self.dr*self.alpha*(self.beta-self.boundary_conds)-self.boundary_conds)+1
        self.bM = self.dt/self.dr*(-(self.sigma**2)/(2*self.dr)-self.alpha*(self.beta-self.boundary_conds))
        self.cM = self.dt*(self.sigma**2)/(6*(self.dr**2))

    @abstractmethod
    def _traverse_grid_(self):
        pass

    @abstractmethod
    def _interpolate_(self):
        return np.interp(self.r0, 
                         self.boundary_conds,
                         self.grid[:, 0])
    
    def price(self):
        self._setup_boundary_conditions_()
        self._setup_coefficients_()
        self._traverse_grid_()
        return self._interpolate_()

"""The Bond finite difference class"""
class FDExplicitBond(FiniteDifferences):
    
    def __init__(self, r0, T, sigma, alpha, beta,rmin, rmax, M, N, couponSchedule, coupon):
        super().__init__(r0, T, sigma, alpha, beta,rmin, rmax, M, N)
        self.couponSchedule = couponSchedule
        self.coupon = coupon

    def _setup_boundary_conditions_(self):
        self.grid[:, -1] = 1 + self.coupon

    def _traverse_grid_(self):
        for j in reversed(self.j_values):
            if j in self.couponSchedule:
                for i in range(self.M+1):
                    if i == 0:                     
                        self.grid[i,j] = self.a0[i]*self.grid[i,j+1] \
                                     +self.b0[i]*self.grid[i+1,j+1] \
                                     +self.c0*self.grid[i+3,j+1] + self.coupon
                    elif i == self.M:
                        self.grid[i,j] = self.aM[i]*self.grid[i,j+1] \
                                     +self.bM[i]*self.grid[i-1,j+1] \
                                     +self.cM*self.grid[i-3,j+1] + self.coupon
                    else:
                        self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] \
                                     +self.b[i]*self.grid[i,j+1] \
                                     +self.c[i]*self.grid[i+1,j+1] + self.coupon
            else:
                for i in range(self.M+1):
                    if i == 0:                     
                        self.grid[i,j] = self.a0[i]*self.grid[i,j+1] \
                                     +self.b0[i]*self.grid[i+1,j+1] \
                                     +self.c0*self.grid[i+3,j+1]
                    elif i == self.M:
                        self.grid[i,j] = self.aM[i]*self.grid[i,j+1] \
                                     +self.bM[i]*self.grid[i-1,j+1] \
                                     +self.cM*self.grid[i-3,j+1]
                    else:
                        self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] \
                                     +self.b[i]*self.grid[i,j+1] \
                                     +self.c[i]*self.grid[i+1,j+1]

"""The TRS finite difference class"""
class FDExplicitTRS(FiniteDifferences):
    
    def __init__(self, r0, T, sigma, alpha, beta,rmin, rmax, M, N, Bond, K, sTRS, is_perfReceiver=True):
        super().__init__(r0, T, sigma, alpha, beta,rmin, rmax, M, N)
        self.Bond = Bond
        self.K = K
        self.sTRS = sTRS
        self.is_perfReceiver = is_perfReceiver
        self.Bond.price()

    def _setup_boundary_conditions_(self):
        index = int(self.Bond.N*self.T/self.Bond.T)
        if self.is_perfReceiver:
            self.grid[:, -1] = self.Bond.grid[:,index] - self.K - self.sTRS*self.T
        else:
            self.grid[:, -1] = -(self.Bond.grid[:,index] - self.K - self.sTRS*self.T)

    def _traverse_grid_(self):
        for j in reversed(self.j_values):
            if j in self.Bond.couponSchedule:
                for i in range(self.M+1):
                    if i == 0:                     
                        self.grid[i,j] = self.a0[i]*self.grid[i,j+1] \
                                     +self.b0[i]*self.grid[i+1,j+1] \
                                     +self.c0*self.grid[i+3,j+1] + self.Bond.coupon
                    elif i == self.M:
                        self.grid[i,j] = self.aM[i]*self.grid[i,j+1] \
                                     +self.bM[i]*self.grid[i-1,j+1] \
                                     +self.cM*self.grid[i-3,j+1] + self.Bond.coupon
                    else:
                        self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] \
                                     +self.b[i]*self.grid[i,j+1] \
                                     +self.c[i]*self.grid[i+1,j+1] + self.Bond.coupon
            else:
                for i in range(self.M+1):
                    if i == 0:                     
                        self.grid[i,j] = self.a0[i]*self.grid[i,j+1] \
                                     +self.b0[i]*self.grid[i+1,j+1] \
                                     +self.c0*self.grid[i+3,j+1]
                    elif i == self.M:
                        self.grid[i,j] = self.aM[i]*self.grid[i,j+1] \
                                     +self.bM[i]*self.grid[i-1,j+1] \
                                     +self.cM*self.grid[i-3,j+1]
                    else:
                        self.grid[i,j] = self.a[i]*self.grid[i-1,j+1] \
                                     +self.b[i]*self.grid[i,j+1] \
                                     +self.c[i]*self.grid[i+1,j+1]

couponSchedule = (25,50,75)
Bond = FDExplicitBond(0.05,2,0.1,0.3,0.01,0.0,0.2,100,100,couponSchedule,0.025)
TRS = FDExplicitTRS(0.05,1,0.1,0.3,0.01,0.0,0.2,100,80,Bond, 1.0,0.02, True)
print(TRS.price())