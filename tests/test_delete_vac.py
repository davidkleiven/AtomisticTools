from ase.build import bulk
from atomtools.ase import delete_vacancies
import unittest
from ase.visualize import view
import random

class TestDeleteVac( unittest.TestCase ):
    def test_del_vac( self ):
        atoms = bulk("Al")*(5,5,5)
        n_vac = 10
        vacs = range(len(atoms))
        random.shuffle(vacs)
        vacs = vacs[:n_vac]
        orig_length = len(atoms)
        for vac in vacs:
            atoms[vac].symbol = "X"

        atoms = delete_vacancies(atoms)
        self.assertEqual( len(atoms), orig_length-n_vac )

        for atom in atoms:
            self.assertNotEqual( atom.symbol, "X" )

if __name__ == "__main__":
    unittest.main()
