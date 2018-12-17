import numpy as np

class DisplacementField(object):
    def __init__(self, ref_cell=None, relaxed_cell=None):
        self.deformation_gradient = self._get_deformation_gradient(ref_cell, relaxed_cell)
        

    def _get_deformation_gradient(self, ref_cell, relaxed_cell):
        F = relaxed_cell.dot(np.linalg.inv(ref_cell))
        return F

    def get(self, x, y, z):
        """Get the displacement position x, y, z."""
        return self.deformation_gradient.dot(np.array([x, y, z]))
