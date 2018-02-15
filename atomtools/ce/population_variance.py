from cemc.wanglandau.ce_calculator import CE
from ase.ce.corrFunc import CorrFunction
import numpy as np
import time
import json

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
        self.calc = CE( self.bc, ecis, initial_cf=self.init_cf )
        self.bc.atoms.set_calculator(self.calc)
        self.elements = self.bc.site_elements[0] # This should be updated to handle different site types
        self.status_every_sec = 30

    def swap_random_atoms( self ):
        indx = np.random.randint(low=0,high=len(self.bc.atoms))
        symb = self.bc.atoms[indx].symbol
        new_symb = symb
        while( new_symb == symb ):
            new_symb = self.elements[np.random.randint(low=0,high=len(self.elements))]
        system_change = [(indx,symb,new_symb)]
        self.bc.atoms._calc.calculate( self.bc.atoms, ["energy"], system_change )

    def estimate( self, n_probe_structures=10000, fname="" ):
        if ( fname != "" ):
            if ( not fname.endswith(".json") ):
                raise ValueError( "The program is going to write a json file so the file extension should be .json" )

        step = 0
        cfs = []
        cur_time = time.time()
        while ( step < n_probe_structures ):
            step += 1
            if ( time.time()-cur_time > self.status_every_sec ):
                print ("Step {} of {}".format(step,n_probe_structures) )
                cur_time = time.time()
            for i in range(len(self.bc.atoms)):
                self.swap_random_atoms()
            new_cfs = self.calc.get_cf()
            cf_array = [new_cfs[key] for key in self.init_cf.keys()] # Make sure that the order is the same as in init_cf
            cfs.append( cf_array )

        cfs = np.array( cfs )
        cov = np.cov( cfs.T )
        mu = np.mean( cfs, axis=0 )

        # Create dictionaries
        keys = self.init_cf.keys()
        mu_dict = {keys[i]:mu[i] for i in range(len(keys))}

        cov_dict = {key:{} for key in keys}
        for i in range(len(keys)):
            for j in range(len(keys)):
                cov_dict[keys[i]][keys[j]] = cov[i,j]

        if ( fname != "" ):
            data = {}
            data["cov"] = cov_dict
            data["mean"] = mu_dict
            with open(fname, 'w') as outfile:
                json.dump( data, outfile )
            print ( "Covariance and mean of the correlation functions are written to {}".format(fname) )
        return cov_dict, mu_dict
