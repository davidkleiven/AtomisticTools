try:
    from ase.ce import evaluate as ev
except ImportError as exc:
    print ( str(exc) )
    print ("Warning! Cannot use the phonon_ce_eval functionality without ASE version that has Cluster Expansion!")
    raise ImportError(exc)

from atomtools.ce import ce_phonon_dos as cpd
import numpy as np
from ase.units import kB
from matplotlib import pyplot as plt

class PhononEval( ev.Evaluate ):
    def __init__(self, BC, cluster_names=None, select_cond=None, lamb=0.0, penalty=None ):
        self.phonon_db = None
        ev.Evaluate.__init__( self, BC, cluster_names=cluster_names, lamb=float(lamb), penalty=penalty )
        self.n_atoms = {}
        self.T = 600

    def _make_cf_matrix( self ):
        """
        Override the make CF matrix. Include only runs that actually
        has a phonon density computed
        """
        manager = cpd.PhononDOS_DB( self.db_name )
        self.phonon_db = manager.get_all()
        atIDs = [res["atID"] for res in self.phonon_db]
        cf_matrix = []
        for i,atID in enumerate(atIDs):
            row = self.db.get( id=atID )
            cf_matrix.append( [row[x]*row.natoms for x in self.cluster_names] )
            dw = self.phonon_db[i]["omega_e"][1] - self.phonon_db[i]["omega_e"][0]
            self.phonon_db[i]["dos_e"] *= 64.0*3.0/(dw*np.sum(self.phonon_db[i]["dos_e"]))
        return np.array(cf_matrix, dtype=float )

    def _get_dft_energy_per_atom(self):
        """
        Override this to return the vibrational free energy per atom
        """
        if ( self.phonon_db is None ):
            manager = cpd.PhononDOS_DB( self.db_name )
            self.phonon_db = manager.get_all()
        energies = []
        for res in self.phonon_db:
            energies.append( self.free_energy( res["omega_e"], res["dos_e"]) )
        self.e_dft = np.array(energies)
        self.e_dft -= self.e_dft[0]
        return True

    def free_energy( self, omega, dos ):
        dw = omega[1]-omega[0]
        beta = 1.0/(kB*self.T)
        return np.sum( dw*np.log(1.0-np.exp(-beta*omega[1:]))*dos[1:])
        return np.sum( dw*np.log(omega[1:])*dos[1:] )