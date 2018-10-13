from ase.clease import Evaluate
import copy
import numpy as np

class EvaluateBootstrap(Evaluate):
    def __init__( self, BC, **kwargs ):
        self.n_bootstrap = kwargs.pop("n_bootstrap")
        super(EvaluateBootstrap,self).__init__(BC,**kwargs)
        self.orig_cfmatrix = copy.deepcopy( self.cf_matrix )
        self._get_dft_energy_per_atom()
        self.orig_e_dft = copy.deepcopy( self.e_dft )

    def create_boot_strap_dataset(self):
        """
        Resamples the dataset
        """
        n_datapoints = len(self.orig_e_dft)
        indx = np.random.randint(low=0,high=n_datapoints,size=n_datapoints)
        for i,indx in enumerate(indx):
            self.e_dft[i] = self.orig_e_dft[indx]
            self.cf_matrix[i,:] = self.orig_cfmatrix[indx,:]

    @property
    def get_eci(self):
        avg_ecis = np.zeros(self.orig_cfmatrix.shape[1])
        for i in range(self.n_bootstrap):
            self.create_boot_strap_dataset()
            super(EvaluateBootstrap,self).get_eci
            avg_ecis += self.eci
        self.eci = avg_ecis/self.n_bootstrap

        # Reset the matrices
        self.cf_matrix = copy.deepcopy(self.orig_cfmatrix)
        self.e_dft = copy.deepcopy(self.orig_e_dft)
        return self.eci

    def _get_eci_loo(self,indx):
        avg_ecis = np.zeros(self.orig_cfmatrix.shape[1])
        for i in range(self.n_bootstrap):
            self.create_boot_strap_dataset()
            eci = super(EvaluateBootstrap,self)._get_eci_loo(indx)
            avg_ecis += eci
        self.cf_matrix = copy.deepcopy( self.orig_cfmatrix )
        self.e_dft = copy.deepcopy( self.e_dft )
        return avg_ecis/self.n_bootstrap
