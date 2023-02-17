from abc import ABC, abstractmethod
import numpy as np

""" 
Base class for sharing 
attributes and functions of FD 
"""
class FiniteDifferences(object):

    def __init__(
        self, S0, K, r=0.05,rTRS=0.05,cyield=0, T=1, 
        Smax=1, M=1, N=1, is_pay_perf=False
    ):
        self.S0 = S0
        self.K = K
        self.r = r
        self.rTRS = rTRS
        self.cyield = cyield
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
    def setup_boundary_conditions(self):
        raise NotImplementedError('Implementation required!')

    @abstractmethod
    def setup_coefficients(self):
        raise NotImplementedError('Implementation required!')

    @abstractmethod
    def traverse_grid(self):
        """  
        Iterate the grid backwards in time 
        """
        raise NotImplementedError('Implementation required!')

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
            self.grid[:,-1] = self.boundary_conds - self.K - self.rTRS*self.K*self.T
        else:
            self.grid[:,-1] = -(self.boundary_conds - self.K - self.rTRS*self.K*self.T)

    def setup_coefficients(self):
        self.a = self.dt*(self.r -self.cyield)*self.i_values
        self.b = 1 - self.dt*self.r

    def traverse_grid(self):
        for j in reversed(self.j_values):
            for i in range(self.M)[2:]:
                self.grid[i,j] = self.a[i]+self.b*self.grid[i,j+1]

TRS = FDExplicitEu(50, 50, r=0.01, rTRS=0.04,cyield=0, T=1., Smax=100, M=100, N=1000, is_pay_perf=True)
print(TRS.price())