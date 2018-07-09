import unittest
from ase.calculators.emt import EMT
from ase.build import bulk
from atomtools.ase import ElasticConstants
from ase.calculators.calculator import Calculator
import numpy as np

# Implement dymmy calculator that supports get_stress
class DummyCalc(Calculator):
    implemented_properties = ["stress"]
    def calculate(self, atoms, properties, system_changes):
        super(DummyCalc, self).calculate(atoms, properties, system_changes)
        self.results = {
            "stress":np.random.rand(6)
        }

class TestElasticConstante(unittest.TestCase):
    def test_no_throw(self):
        no_throw = True
        msg = ""
        try:
            atoms = bulk("Al")
            calc = DummyCalc()
            atoms.set_calculator(calc)
            el = ElasticConstants(atoms)
            C = el.get()
            B = el.bulk_modulus()
            G = el.shear_modulus()
            poisson = el.poisson_ratio()
        except Exception as exc:
           msg = str(exc)
           no_throw = False
        self.assertTrue(no_throw, msg=msg)

if __name__ == "__main__":
    unittest.main()
