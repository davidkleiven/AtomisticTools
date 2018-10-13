from ase.clease import Evaluate

class EvaluateGaussianEOFPrior(Evaluate):
    def __init__( self, **kwargs ):
        if ( "prior_weights" not in kwargs.keys() ):
            raise ValueError( "Prior weights has to be given" )
        self.prior_lambda = kwargs.pop("prior_weights")
        Evaluate.__init__(self, **kwargs)
        self._current_lambda_indx = 0

    def get_current_weight(self):
        try:
            # User supplied an array. Return the one of them
            return self.prior_lambda[self._current_lambda_indx]
        except:
            # User supplied only one value for lambda. Return it.
            return self.prior_lambda

    def _get_dft_energy_per_atom(self):
        """
        Override in parents class
        """
        e_dft = []
        for row in self.setting.db.select(self.select_cond):
            is_gaussian_eof_prior = row.get("gaussian_eof_prior",default=False)
            if ( is_gaussian_eof_prior ):
                ref_en_mix_new = row.mix_energy_new
                ref_en_mix_old = row.mix_energy_old
                E = row.energy - ( ref_en_mix_old-ref_en_mix_new )
                e_dft.append( E*np.sqrt(self.get_current_weight()) )
            else:
                e_dft.append( row.energy/row.natoms )
        self.e_dft = np.array(e_dft)
