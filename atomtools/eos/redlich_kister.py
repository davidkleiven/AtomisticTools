import numpy as np

class RedlichKister( object ):
    def __init__( self, ref_energies ):
        self.ref_energies = ref_energies
        self.coefficients = None
        self.order = None

    def _check_composition( self, comp ):
        """
        Checks that the composition argument is OK
        """
        if ( len(composition.keys()) > 2 ):
            raise NotImplementedError( "Currently only binary alloys are supported" )

        for key in comp.keys():
            if ( key not in self.ref_energies.keys() ):
                msg = "The composition argument has to have the following keys"
                msg += "{}. Given: {}".format(self.ref_energies.keys(), comp.keys())
                raise ValueError( msg )

        ref_length = comp[comp.keys()[0]]
        for key,value in comp.iteritems():
            if ( len(value) != ref_length ):
                raise ValueError( "All elements has to have the same number of datapoints" )

    def free_energy2excess_energy( self, free_energy, composition ):
        """
        Converts the Free Energy to Excess energy
        """
        for key,value in composition.iteritems():
            free_energy -= self.ref_energies[key]*value
        return free_energy

    def fit( self, free_energy, composition, order=3 ):
        """
        Perform fit
        """
        n_data_points = len(composition[composition.keys()[0]])
        A = np.zeros((n_data_points,order))
        xa = np.array( composition[composition.keys()[0]] )
        xb = np.array( composition[composition.keys()[1]] )

        for power in range(order):
            A[:,power] = (xa-xb)**power
        A *= xa*xb
        free_energy = self.free_energy2excess_energy(free_energy,composition)
        self.coefficients,rank,s = np.linalg.lstsq(A,free_energy)
        self.order = order
        return self.coefficients

    def eval( self, composition ):
        """
        Evaluates the enthalpy of formation
        """
        if ( self.coefficients is None ):
            msg = "The fitting coefficients have not been computed."
            msg += " These can be found by invoking the fit function"
            raise ValueError( msg )

        excess = 0.0
        xa = composition[composition.keys()[0]]
        xb = composition[composition.keys()[1]]
        for power in self.order:
            excess += (xa-xb)**power
        excess *= xa*xb
        return excess
