import unittest
import numpy as np
try:
    from atomtools.ase import DisplacementField
    available = True
    reason = ""
except ImportError as exc:
    available = False
    reason = str(exc)


class TestDisplacementField(unittest.TestCase):
    def test_corner_disp(self):
        cell1 = np.array([(0.1, 0, 0.1), (0, 0.2, 0.1), (0.2, 0.3, 0)])
        F = np.array([(1.02, 0.001, 0.002), (0.001, 1.02, 0.003), (0.002, 0.003, 1.03)])

        # Make sure strain is symmetric
        cell2 = F.dot(cell1)

        disp_field = DisplacementField(ref_cell=cell1, relaxed_cell=cell2)

        # Check that the diagonal is correct
        diag1 = np.sum(cell1, axis=1)
        diag2 = np.sum(cell2, axis=1)
        displacement = diag2 - diag1
        disp_calc = disp_field.get(diag1[0], diag1[1], diag1[2]) - diag1
        
        # Correct by construction
        self.assertTrue(np.allclose(disp_calc, displacement))

        # Check the other corners
        disp1 = cell2[:, 0] - cell1[:, 0]
        disp_calc = disp_field.get(cell1[0, 0], cell1[1, 0], cell1[2, 0])
        self.assertTrue(np.allclose(cell2[:, 0], disp_calc))

if __name__ == "__main__":
    unittest.main()

