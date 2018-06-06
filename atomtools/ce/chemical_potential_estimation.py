import numpy as np
from itertools import combinations
from scipy.misc import comb
from matplotlib import pyplot as plt

class ChemicalPotentialEstimator(object):
    def __init__( self, singlets=None, energies=None ):
        self.singlets = singlets
        self.energies = energies
        self.coeff = None
        self.term_list = []
        self.term_type = []

    def eval( self, x ):
        if ( self.coeff is None ):
            self.parabolic_fit()

        dim = self.singlets.shape[1]
        value = self.coeff[0]
        n_linear_terms = dim
        n_quadratic_terms = dim
        n_bilinear_terms = comb(dim,2)
        counter = 1
        for i in range(0,n_linear_terms):
            value += self.coeff[counter]*x[i]
            counter += 1
        for indx in combinations(range(dim),2):
            value += self.coeff[counter]*x[indx[0]]*x[indx[1]]
            counter += 1
        for i in range(n_quadratic_terms):
            value += self.coeff[counter]*x[i]**2
            counter += 1
        return value

    def parabolic_fit( self ):
        """
        Perform a parabolic fit
        """
        dim = self.singlets.shape[1]
        n_linear_terms = dim
        n_quadratic_terms = dim
        n_bilinear_terms = comb(dim,2)
        n_terms = 1 + n_linear_terms + n_bilinear_terms + n_quadratic_terms
        X = np.ones( (len(self.energies), int(n_terms) ) )
        counter = 1
        self.term_list = []
        self.term_list.append( (None,) )
        self.term_type = []
        self.term_type.append( "constant" )
        for i in range(0,n_linear_terms):
            X[:,counter] = self.singlets[:,i]
            counter += 1
            self.term_list.append( (i,) )
            self.term_type.append( "linear" )

        for indx in combinations(range(dim),2):
            X[:,counter] = self.singlets[:,indx[0]]*self.singlets[:,indx[1]]
            counter += 1
            self.term_list.append( (indx[0],indx[1]) )
            self.term_type.append( "bilinear" )

        for i in range(n_quadratic_terms):
            X[:,counter] = self.singlets[:,i]**2
            counter += 1
            self.term_list.append( (i,) )
            self.term_type.append( "quadratic" )

        self.coeff, res, rank, s = np.linalg.lstsq( X, self.energies )

    def plot(self):
        """
        Creates a plot of the qualitty of the fit
        """
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.energies, self.energies )
        fitted = []
        for i in range(self.singlets.shape[0]):
            fitted.append( self.eval(self.singlets[i,:]) )
        ax.plot( self.energies, fitted, "o", mfc="none" )

    def deriv( self, x, direction ):
        """
        Computes the derivative in certain direction
        """
        if ( self.coeff is None ):
            self.parabolic_fit()

        value = 0.0
        for i in range(len(self.term_list)):
            if ( direction in self.term_list[i] ):
                if ( self.term_type[i] == "linear" ):
                    value += self.coeff[i]
                elif ( self.term_type[i] == "quadratic" ):
                    value += 2.0*self.coeff[i]*x[direction]
                elif ( self.term_type[i] == "bilinear" ):
                    # Get the other index
                    indx = direction
                    for entry in self.term_list[i]:
                        if ( entry != direction ):
                            indx = entry
                            break
                    value += self.coeff[i]*x[indx]
        return value
