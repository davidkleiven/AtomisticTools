from matplotlib import pyplot as plt

class EquationOfState(object):
    def __init__( self, volume, energy ):
        self.volume = volume
        self.energy = energy

    def evaluate( self, V ):
        raise NotImplementedError( "This function has to be implemented in subclasses" )

    def deriv( self, V ):
        raise NotImplementedError( "Derivatives has to be implemented in subclasses" )

    def double_deriv( self, V ):
        raise NotImplementedError( "Double derivative has to be implemented in subclasses" )

    def plot( self, latex=False ):
        """
        Plots the result
        """
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.volume, self.energy, "o", mfc="none" )
        V = np.linspace( 0.95*np.min(self.volume), 1.05*np.max(self.volume), 100 )
        fitted = self.evaluate(V)
        ax.plot( V, fitted )
        if ( latex ):
            ax.set_xlabel( "Volume (\$\SI{}{\\angstrom^3}\$)")
        else:
            ax.set_xlabel( r"Volume (\AA^3)" )
        ax.set_ylabel( "Energy (eV)" )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        return fig

    def bulk_modulus( self, V ):
        """
        Computes the bulk modulus as a function of volume
        """
        B = -self.double_deriv(V)/V
        return B
