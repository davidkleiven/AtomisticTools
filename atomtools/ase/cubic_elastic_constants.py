from ase.db import connect
from ase.constaints import StrainFilter
from ase.optimize import BFGS

class ElasticConstants(object):
    def __init__( self, db_name=None, atoms=None ):
        self.orig_atoms = self.atoms.copy()
        self.orig_atoms.set_calculator( self.atoms._calc )

    def run_relaxations( self, vol_scale_factor=1.05, smax=0.003 ):
        """
        Runs 6 relaxations
        """
        strain_components = ["xx","yy","zz","yz","xz","xy"]
        for i in range(6):
            atoms = self.orig_atoms.copy()
            atoms.set_calculator( self.orig_atoms._calc )
            mask = np.zeros(6)
            mask[i] = 1
            strfilter = StrainFilter( atoms, mask=mask )
            relaxer = BFGS(strfilter)
            relaxer.run( fmax=smax*atoms.get_volume() )
            db = connect( self.db_name )
            db.write( atoms, key_value_pairs={"non_zero_strain":strain_components[i]} )

    def fit_elastic_constants( self ):
        """
        Fit the elastic constants
        """
        db = connect( self.db_name )
        stresses = []
        strains = []
        for row in db.select():
            atoms = row.toatoms()
            try:
                stress = atoms.get_stress()
                strain = atoms.get_strain()
            except Exception as exc:
                print (str(exc))
                continue
            stresses.append(stress)
            strains.append(strain)
