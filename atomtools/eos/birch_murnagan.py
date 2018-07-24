from atomtools.eos.equation_of_state import EquationOfState
import numpy as np

class BirschMurnagan( EquationOfState ):
    def __init__( self, volume, energy ):
        super(BirschMurnagan,self).__init__(volume,energy)
        self.a = None
        self.b = None
        self.c = None
        self.d = None

    def fit( self ):
        """
        Fits the 4 parameters in the equation of state to the
        equation of state
        """
        x = np.zeros((len(self.volume),5))
        x[:,0] = 1.0
        x[:,1] = self.volume**(-1.0/3.0)
        x[:,2] = self.volume**(-2.0/3.0)
        x[:,3] = self.volume**(-1.0)
        res, residuals,rank,s = np.linalg.lstsq( x, self.energy )
        self.a = res[0]
        self.b = res[1]
        self.c = res[2]
        self.d = res[3]

    def perform_fit(self):
        return self.a is None or self.b is None or self.c is None or self.d is None

    def evaluate( self, volume ):
        """
        Evaluates the energy using the fitted parameters
        """
        if ( self.perform_fit() ):
            self.fit()
        return self.a + self.b*volume**(-1.0/3.0) + self.c*volume**(-2.0/3.0) + self.d/volume

    def deriv( self, V ):
        """
        Evaluates the derivative with respect to volume
        """
        if ( self.perform_fit() ):
            self.fit()
        return -(self.b/3.0)*V**(-4.0/3.0) - (2.0*self.c/3.0)*V**(-5.0/3.0) - d*V**(-2.0)

    def double_deriv( self, V ):
        """
        Evaluates the double derivative with respect to volume
        """
        if ( self.perform_fit() ):
            self.fit()
        return (4.0*self.b/9.0)*V**(-7.0/3.0) + (10.0*self.c/9.0)*V**(-8.0/3.0) + 2.0*self.d*V**(-3.0)
