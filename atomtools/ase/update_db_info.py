
from ase.db import connect

class UpdateDBInfo(object):
    """
    Updates information in a database of atoms
    """
    def __init__( self, db_name, uid, atoms ):
        """
        Updates the database with info during run
        """
        self.db_name = db_name
        self.uid = uid
        self.atoms = atoms

    def __call__( self ):
        energy = self.atoms.get_potential_energy()
        stress = 0.0
        try:
            stresses = self.atoms.get_stress()
            stress = np.max(stresses)
        except:
            pass

        fmax = 0.0
        try:
            force = self.atoms.get_forces()
            f = np.sqrt( np.sum( force**2, axis=1 ) )
            fmax = np.max(f)
        except:
            pass

        db = connect( self.db_name )
        row = db.get( id=self.uid )
        kvp = row.key_value_pairs
        update = True
        if ( "min_energy" in kvp.keys() ):
            update = ( energy < kvp["min_energy"] )

        if ( update ):
            db.update( self.uid, min_energy=energy, max_stress=stress, max_force=fmax )
