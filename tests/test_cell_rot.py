import unittest
from atomtools.ase import rotate_atoms_and_cell
from ase import build

class TestRotateAtoms( unittest.TestCase ):
    def test_rotate(self):
        atoms = build.bulk("Al")
        target_z_dirs = [(0,0,1),(0,1,0),(1,0,0),(3,1,1),(-2,1,3),(1,1,1),(1,1,2),(7,1,1)]
        msg = ""
        no_throw = True
        for zdir in target_z_dirs:
            try:
                rotate_atoms_and_cell( atoms, target_z_direction=zdir )
            except RuntimeError as exc:
                msg = str(exc)
                no_throw = False
                break

        self.assertTrue( no_throw, msg=msg )

if __name__ == "__main__":
    unittest.main()
