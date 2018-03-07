try:
    from ase.ce import evaluate as ev
    from ase.ce.corrFunc import CorrFunction
except ImportError as exc:
    print ( str(exc) )
    print ("Warning! Cannot use the phonon_ce_eval functionality without ASE version that has Cluster Expansion!")
    raise ImportError(exc)

from atomtools.ce import ce_phonon_dos as cpd
import numpy as np
from ase.units import kB
from matplotlib import pyplot as plt
from scipy.integrate import simps
from atomtools.eos.birch_murnagan import BirschMurnagan
from ase.db import connect
from matplotlib import pyplot as plt

class PhononEvalEOS( ev.Evaluate ):
    def __init__(self, BC, phonon_dbname, cluster_names=None, select_cond=None, lamb=0.0, penalty=None ):
        self.ph_db = connect( phonon_dbname )
        self._temperature = 600
        self.volume = {}
        self.energy = {}
        self.atoms_count = {}
        self.tot_number_of_atoms = {}
        self.BC = BC
        ev.Evaluate.__init__( self, BC, cluster_names=cluster_names, lamb=float(lamb), penalty=penalty )

    def _make_cf_matrix( self ):
        """
        Override the CF matrix
        """
        cf_matrix = []
        cf = CorrFunction( self.BC )
        for row in self.ph_db.select():
            gid = row.key_value_pairs["groupID"]
            if ( gid in self.volume.keys() ):
                self.volume[gid].append( row.volume )
                self.energy[gid].append( row.energy )
            else:
                self.volume[gid] = [row.volume]
                self.energy[gid] = [row.energy]
                self.atoms_count[gid] = row.count_atoms()
                self.tot_number_of_atoms[gid] = row.natoms
                atoms = self.ph_db.get_atoms(id=row.id)
                if ( len(atoms) == 1 ):
                    atoms = atoms*(4,4,4)
                elif ( len(atoms) == 8 ):
                    atoms = atoms*(2,2,2)
                else:
                    raise ValueError( "Can only handle cases where the atoms object contains 1 or 8 atoms!" )
                atoms.set_cell( self.BC.atoms.get_cell(), scale_atoms=True )
                corr_funcs = cf.get_cf_by_cluster_names( atoms, self.cluster_names )
                cf_matrix.append( [corr_funcs[x] for x in self.cluster_names] )
        return np.array( cf_matrix, dtype=float )

    def _get_dft_energy_per_atom(self):
        """
        Override the function in the main class to return the Free Energy associated
        to elastic and vibrational energy
        """
        e_dft = []
        for key in self.atoms_count.keys():
            eos = BirschMurnagan( np.array(self.volume[key]), np.array(self.energy[key]) )
            eos.set_average_mass( self.atoms_count[key] )
            fvib_el = eos.beta_elastic_vib_free_energy( [self.temperature], natoms=self.tot_number_of_atoms[key] )
            e_dft.append(fvib_el)
        self.e_dft = np.array( e_dft )
        self.e_dft -= np.min(self.e_dft)
        return True

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self,T):
        self.e_dft = None # Set the e_dft to None as they have to be computed again when the temperature change
        self.eci = None
        self._temperature = T

class PhononEval( ev.Evaluate ):
    def __init__(self, BC, cluster_names=None, select_cond=None, lamb=0.0, penalty=None, filters=[], hte=True ):
        self.phonon_db = None
        self.formulas = []
        self.filters = filters
        self.counts = []
        ev.Evaluate.__init__( self, BC, cluster_names=cluster_names, lamb=float(lamb), penalty=penalty )
        self.n_atoms = {}
        self.T = 600
        self.hte = hte

    def include_structure( self, count ):
        for f in self.filters:
            if ( not f(count) ):
                return False
        return True

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
            self.counts.append( row.count_atoms() )
            if ( self.include_structure(self.counts[-1]) ):
                self.formulas.append(row.formula)
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
        for i,res in enumerate(self.phonon_db):
            if ( self.include_structure(self.counts[i]) ):
                if ( self.hte ):
                    energies.append( self.logw(res["omega_e"], res["dos_e"]) )
                else:
                    energies.append( self.free_energy( res["omega_e"], res["dos_e"]) )
        self.e_dft = np.array(energies)
        print (self.formulas)
        print (self.e_dft)
        #self.e_dft -= self.e_dft[0]
        return True

    def logw( self, omega, dos ):
        dw = omega[1]-omega[0]
        return simps( dw*np.log(omega[1:])*dos[1:] )

    def free_energy( self, omega, dos ):
        dw = omega[1]-omega[0]
        beta = 1.0/(kB*self.T)
        return simps( dw*np.log(1.0-np.exp(-beta*omega[1:]))*dos[1:])
        return np.sum( dw*np.log(omega[1:])*dos[1:] )
