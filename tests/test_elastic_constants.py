"""Unit tests for the elastic constants."""
import unittest
from ase.build import bulk
from atomtools.ase import ElasticConstants
from ase.calculators.calculator import Calculator
import numpy as np
import os


class DummyCalc(Calculator):
    """Dummy calculator for extracting stress."""

    implemented_properties = ["stress"]

    def calculate(self, atoms, properties, system_changes):
        """Override calculator method."""
        Calculator.calculate(self, atoms, properties, system_changes)
        self.results = {
            "stress": np.random.rand(6)
        }


db_name = "test_elastic.db"


class TestElasticConstante(unittest.TestCase):
    """Unit test case for elastic constants."""

    def test_no_throw(self):
        """Test that code runs without throwing exceptions."""
        no_throw = True
        msg = ""
        try:
            atoms = bulk("Al")
            calc = DummyCalc()
            atoms.set_calculator(calc)
            el = ElasticConstants(atoms, db_name)
            el.prepare_db()
            for uid in range(1, 25):
                el.run(uid, calc)
            el.get()
            el.bulk_modulus()
            el.shear_modulus()
            el.poisson_ratio
            cell = atoms.get_cell()
            ElasticConstants.get_strain(cell, cell, principal=True)
            el.isotropic_elastic_tensor
            os.remove(db_name)
        except Exception as exc:
            msg = str(exc)
            no_throw = False
        self.assertTrue(no_throw, msg=msg)

    def test_mandel_conversion_rank4(self):
        atoms = bulk("Al")
        calc = DummyCalc()
        atoms.set_calculator(calc)
        el = ElasticConstants(atoms, db_name)

        tensor = np.random.rand(6, 6)
        full = el._to_full_rank4(tensor)
        tensor_orig = el._to_mandel_rank4(full)
        self.assertTrue(np.allclose(tensor, tensor_orig))

if __name__ == "__main__":
    unittest.main()
