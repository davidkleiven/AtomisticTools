import unittest
from ase.build import bulk
from atomtools.ce import BondLengthDistribution

class TestBondLength(unittest.TestCase):
    def test_no_throw(self):
        from random import choice
        atoms = bulk("Al")*(4, 4, 4)
        elements = ["Mg", "Si", "Al"]
        for atom in atoms:
            atom.symbol = choice(elements)
        bond_lengths = BondLengthDistribution(atoms, num_bonds=12)
        blengths = bond_lengths.calculate_bond_lengths()
        bond_lengths.show_bond_lengths(show=False)
        keys = list(blengths.keys())
        expected_keys = sorted(["Al-Al", "Al-Mg", "Al-Si", "Mg-Si",
                                "Mg-Mg", "Si-Si"])
        self.assertTrue(sorted(keys) == expected_keys)

if __name__ == "__main__":
    unittest.main()
