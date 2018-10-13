from ase.clease.evaluate import Evaluate
import numpy as np

class EvaluateDeviation( Evaluate ):
    def __init__( self, *args, **kwargs ):
        Evaluate.__init__( self, *args, **kwargs )
        self.ref_energies = {}
        self.model_energies = []

    def locate_reference_energies(self):
        for row in self.db.select(self.select_cond):
            count = row.count_atoms()
            if ( len(count.keys()) == 1 ):
                self.ref_energies[count.keys()[0]] = row.energy/row.natoms

    def _get_dft_energy_per_atom(self):
        self.locate_reference_energies()
        e_dft = []
        for row in self.db.select(self.select_cond):
            count = row.count_atoms()
            deviation = row.energy
            mod_energy = 0.0
            for key,value in count.iteritems():
                try:
                    mod_energy += self.ref_energies[key]*value
                except KeyError as exc:
                    msg = "The element does not have any reference energy in the database."
                    msg += "It does not make sense to use EvaluateDeviation without reference energies (energy pure elements)"
                    msg += "Use Evaluate instead"
                    raise KeyError(msg)

            deviation -= mod_energy
            self.model_energies.append(mod_energy/row.natoms)
            e_dft.append( deviation/row.natoms )
        self.e_dft = np.array(e_dft)
        return True
