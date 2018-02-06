from wanglandau.ce_calculator import CE
from ase.ce.corrFunc import CorrFunction
import numpy as np

class PopulationVariance( object ):
    """
    Object that estimates the covariance and the mean of all the
    correlation functions in the population
    """
    def __init__( self, BC ):
        self.bc = BC
        cf_obj = CorrFunction(self.bc)
        self.init_cf = cf_obj.get_cf( self.bc.atoms )
        ecis = {key:1.0 for key in self.init_cf.keys()} # ECIs do not matter here
        self.calc = CE( self.bc, ecis, initial_cf=init_cf )
        self.bc.atoms.set_calculator(self.calc)
        self.elements = self.bc.site_elements[0] # This should be updated to handle different site types

    def swap_random_atoms( self ):
        indx = np.random.randint(low=0,high=len(self.bc.atoms))
        symb = self.bc.atoms[indx].symbol
        new_symb = symb
        while( new_symb == symb ):
            new_symb = self.elements[np.random.randint(low=0,high=len(self.elements))]
        system_change = [(indx,symb,new_symb)]
        self.bc.atoms._calc.calculate( self.bc.atoms, ["energy"], system_change )

    def estimate( self, n_probe_structures=10000 ):
        step = 0
        cfs = []
        while ( step < n_probe_structures ):
            step += 1
            for i in range(len(self.bc.atoms)):
                self.swap_random_atoms()
                cfs.append( self.calc.get_cf() )

        cfs = np.array( cfs )
        cov = np.cov( cfs )
        mu = np.mean( cfs, axis=0 )
        return cov, mu
