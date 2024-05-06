# 
# lattice.py                                                               
# 
# D. Clarke
#
# A basic lattice class for use with statistical mechanical models. 
# 


import numpy as np
import itertools
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType


class Lattice():


    def __init__(self, geometry, example):

        checkType(geometry,"array")

        if np.asarray(example).shape:
            self.nodeShape = np.asarray(example).shape
            self.grid = np.zeros((*geometry, *self.nodeShape))
        else:
            self.nodeShape = 1 
            self.grid = np.zeros(geometry)

        self.Nd       = len(geometry)
        self.geometry = np.array(geometry)
        self.vol      = np.prod(self.geometry)
        self.ranges   = [np.array(range(num)) for num in self.geometry]
        self.bulk     = list(itertools.product(*self.ranges))

        logger.details(f'Initialized {self.Nd}d lattice.')
        logger.details('        dims = ', geometry)
        logger.details('  grid shape = ', self.grid.shape)


    def march(self,coord,mu,step) -> np.ndarray:
        """ Move forward from coordinate coord in direction mu a distance step.

        Args:
            coord (np.array)
            mu (int)
            step (int)

        Returns:
            np.array: New coordinate. 
        """
        newcoord = np.array(coord)
        np.add.at(newcoord, mu, step)
        return newcoord % self.geometry


    def setElement(self,coord,RHS):
        self.grid[tuple(coord)]=RHS


    def getElement(self,coord):
        return self.grid[tuple(coord)]


    # define iterateOverDomain('bulk','random',etc)
    def iterateOverBulk(self,func) -> np.ndarray:
        """ Carry out function func on every site of the lattice. 

        Args:
            func (func): Action to be carried out on every site. Must take coord as
                         its only argument.
        """
        container=[]
        for coord in self.bulk:
            container.append(func(coord))
        return np.array(container)


    def bulkReduce(self,func) -> float:
        return np.sum(self.iterateOverBulk(func))/self.vol


    def iterateOverRandom(self,func) -> np.ndarray:
        """ Carry out function func on every site of the lattice. 

        Args:
            func (func): Action to be carried out on every site. Must take coord as
                         its only argument.
        """
        permutation = np.random.choice(range(self.vol),size=self.vol,replace=False)
        container=[]
        for i in permutation:
            func(self.bulk[i])
        return np.array(container)

