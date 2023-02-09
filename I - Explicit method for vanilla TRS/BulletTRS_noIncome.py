from abc import ABC, abstractmethod
import numpy as np

""" 
Base class for sharing 
attributes and functions of FD 
"""
class FiniteDifferences(object):

    def __init__(
        self, S0, K, r=0.05,r0=0.05, T=1, 
        Smax=1, M=1, N=1, is_pay_perf=False
    ):
        self.S0 = S0
        self.K = K
        self.r = r
        self.r0 = r0
        self.T = T
        self.Smax = Smax
        self.M, self.N = M, N
        self.is_receive_perf = not is_pay_perf

        self.i_values = np.arange(self.M)
        self.j_values = np.arange(self.N)
        self.grid = np.zeros(shape=(self.M+1, self.N+1))
        self.boundary_conds = np.linspace(0, Smax, self.M+1)

    @property
    def dS(self):
        return self.Smax/float(self.M)

    @property
    def dt(self):
        return self.T/float(self.N)

    @abstractmethod
    def interpolate(self):
        """
        Use piecewise linear interpolation on the initial
        grid column to get the closest price at S0.
        """
        return np.interp(
            self.S0, self.boundary_conds, self.grid[:,0])

    def price(self):
        self.setup_boundary_conditions()
        self.setup_coefficients()
        self.traverse_grid()
        return self.interpolate()

""" 
Explicit method of Finite Differences 
"""
class FDExplicitEu(FiniteDifferences):

    def setup_boundary_conditions(self):
        if self.is_receive_perf:
            self.grid[:,-1] = self.boundary_conds - self.K - self.r0*self.K*self.T
            self.grid[-1,:-1] = (self.Smax-self.K - self.r0*self.K*self.T) * \
                np.exp(-self.r*self.dt*(self.N-self.j_values))
        else:
            self.grid[:,-1] = -(self.boundary_conds - self.K - self.r0*self.K*self.T)
            self.grid[0,:-1] = (self.K-self.Smax + self.r0*self.K*self.T) * \
                np.exp(-self.r*self.dt*(self.N-self.j_values))

    def setup_coefficients(self):
        self.a = 0.5*self.dt*self.r*self.i_values
        self.b = 1 - self.dt*self.r
        self.c = 0.5*self.dt*self.r*self.i_values

    def traverse_grid(self):
        for j in reversed(self.j_values):
            for i in range(self.M)[2:]:
                self.grid[i,j] = \
                    self.a[i]*self.grid[i-1,j+1] +\
                    self.b*self.grid[i,j+1] + \
                    self.c[i]*self.grid[i+1,j+1]

TRS = FDExplicitEu(50, 50, r=0.1, r0=0.04, T=5./12., Smax=100, M=100, N=1000, is_pay_perf=True)
print(TRS.price())