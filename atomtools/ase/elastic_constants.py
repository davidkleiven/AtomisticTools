"""Module for calculating elastic constants."""
from ase.db import connect
import numpy as np


class ElasticConstants(object):
    """Class that estimate the elastic parameters."""

    def __init__(self, atoms, db_name):
        """Class esitimating elastic constants.

        :param atoms: Relaxed atoms object
        :param db_name: Database name of structures that needs to be computed
        """
        self.atoms = atoms

        self.delta_no_shear = [-0.01, -0.005, 0.005, 0.01]
        self.delta_shear = [-0.06, -0.03, 0.03, 0.06]
        self.data = []
        self.elastic_tensor = None
        self.db_name = db_name

    def _to_voigt(self, tensor):
        """Convert tensor to Voigt notation."""
        voigt = np.zeros(6)
        voigt[0] = tensor[0, 0]
        voigt[1] = tensor[1, 1]
        voigt[2] = tensor[2, 2]
        voigt[3] = tensor[1, 2]
        voigt[4] = tensor[0, 2]
        voigt[5] = tensor[0, 1]
        return voigt

    def _compute_no_shear(self):
        """Compute the energy for the non shear configurations."""
        identity = np.identity(3)
        db = connect(self.db_name)
        strain_type = 0
        for delta in self.delta_no_shear:
            for i in range(3):
                atoms = self.atoms.copy()
                atoms.set_calculator(self.atoms.get_calculator())
                F = np.identity(3)
                F[i, i] = 1 + delta
                strain = 0.5*(F.T.dot(F) - identity)
                cell = atoms.get_cell()  # NOTE: not transpose by purpose

                # Scale i-th component of each lattice vector
                cell = cell.dot(F)
                atoms.set_cell(cell, scale_atoms=True)
                kvp = {"strain_type": strain_type}
                strain = self._to_voigt(strain)
                db.write(atoms, data={"strain": strain}, key_value_pairs=kvp)

    def _compute_shear(self):
        """Compute the stresses for sheared configurations."""
        identity = np.identity(3)
        element = [(0, 1), (0, 2), (1, 2)]
        db = connect(self.db_name)
        strain_type = 3
        for delta in self.delta_shear:
            for e in element:
                atoms = self.atoms.copy()
                atoms.set_calculator(self.atoms.get_calculator())
                F = np.identity(3)
                F[e[0], e[1]] = delta
                strain = 0.5*(F.T.dot(F) - identity)
                cell = atoms.get_cell()
                cell = cell.dot(F)
                atoms.set_cell(cell, scale_atoms=True)
                kvp = {"strain_type": strain_type}
                strain = self._to_voigt(strain)
                db.write(atoms, data={"strain": strain}, key_value_pairs=kvp)
                strain_type += 1

    def prepare_db(self):
        """Prepare database for DFT calculations.

        Puts entries into the database which needs to be evaluated with
        DFT
        """
        self._compute_no_shear()
        self._compute_shear()

    def run(self, uid, calc):
        """Run one job."""
        db = connect(self.db_name)
        row = db.get(id=uid)
        strain = row.data["strain"]
        kvp = row.key_value_pairs
        atoms = db.get_atoms(id=uid)
        atoms.set_calculator(calc)
        stress = atoms.get_stress()
        del db[uid]
        db.write(atoms, data={"stress": stress, "strain": strain},
                 key_value_pairs=kvp)

    def get(self, select_cond=[]):
        """Compute the elastic properties."""
        db = connect(self.db_name)
        self.data = []
        for row in db.select():
            d = row.data
            d["strain_type"] = row.key_value_pairs["strain_type"]
            self.data.append(d)

        stress_matrix = np.zeros((6, len(self.data)))
        strain_matrix = np.zeros_like(stress_matrix)
        for i, d in enumerate(self.data):
            stress_matrix[:, i] = d["stress"]
            strain_matrix[:, i] = d["strain"]

        self.elastic_tensor = stress_matrix.dot(np.linalg.pinv(strain_matrix))
        return self.elastic_tensor

    @property
    def compliance_tensor(self):
        """Return the compliance tensor."""
        return np.linalg.inv(self.elastic_tensor)

    def _bulk_mod_voigt(self):
        """Bulk modulus by the Voigt average."""
        C = self.elastic_tensor
        Kv = C[0, 0] + C[1, 1] + C[2, 2]
        Kv += 2.0*(C[0, 1] + C[1, 2] + C[2, 0])
        Kv /= 9.0
        return Kv

    def _bulk_mode_reuss(self):
        """Bulk modulus by the reuss average."""
        S = self.compliance_tensor
        inv_KR = S[0, 0] + S[1, 1] + S[2, 2]
        inv_KR += 2.0*(S[0, 1] + S[1, 2] + S[2, 0])
        kR = 1.0/inv_KR
        return kR

    def _shear_mod_voigt(self):
        """Shear modulus by Voigt average."""
        C = self.elastic_tensor
        Gv = C[0, 0] + C[1, 1] + C[2, 2]
        Gv -= (C[0, 1] + C[1, 2] + C[2, 0])
        Gv += 3.0*(C[3, 3] + C[4, 4] + C[5, 5])
        return Gv/15.0

    def _shear_mod_reuss(self):
        """Shear modulus by Reuss average."""
        S = self.compliance_tensor
        inv_Gr = 4.0*(S[0, 0] + S[1, 1] + S[2, 2])
        inv_Gr -= 4.0*(S[0, 1] + S[1, 2] + S[2, 0])
        inv_Gr += 3.0*(S[3, 3] + S[4, 4] + S[5, 5])
        Gr = 1.0/inv_Gr
        return Gr/15.0

    def shear_modulus(self, mode="VRH"):
        """Compute the shear modulus."""
        allowed_modes = ["VRH", "V", "R"]
        if mode not in allowed_modes:
            raise ValueError("Mode has to be one of {}".format(allowed_modes))

        if mode == "V":
            return self._shear_mod_voigt()
        if mode == "R":
            return self._shear_mod_reuss()

        Gv = self._shear_mod_voigt()
        Gr = self._shear_mod_reuss()
        return 0.5*(Gv + Gr)

    def bulk_modulus(self, mode="VRH"):
        """Compute the bulk modulus."""
        allowed_modes = ["VRH", "V", "R"]
        if mode not in allowed_modes:
            raise ValueError("Mode has to be one of {}".format(allowed_modes))
        if mode == "V":
            return self._bulk_mod_voigt()
        elif mode == "R":
            return self._bulk_mode_reuss()

        Kv = self._bulk_mod_voigt()
        Kr = self._bulk_mode_reuss()
        return 0.5*(Kv + Kr)

    def youngs_modulus(self, mode="VRH"):
        """Compute Young's modulus."""
        B = self.bulk_modulus(mode=mode)
        return 3.0*B*(1.0 - 2.0*self.poisson_ratio)

    @property
    def poisson_ratio(self):
        """Compute the isotropic Poisson ratio."""
        K_vrh = self.bulk_modulus(mode="VRH")
        G_vrh = self.shear_modulus(mode="VRH")
        return (3.0*K_vrh - 2.0*G_vrh)/(6.0*K_vrh + 2.0*G_vrh)

    @staticmethod
    def get_strain(ref_cell, strained_cell, principal=False):
        """Compute the eigenstrain of the strained cell.

        The eigenstrain is defined as the strain of strained_cell if put
        into the shape of ref_cell. Hence,
        ref_cell = (I + <strain tensor>) * <strained_cell>

        NOTE: Cell vectors are assumed to by column vectors of the two cells
        :param ref_cell: Reference cell 3x3 matrix
        :param strained_cell: Cell to be deformed 3x3 matrix
        """
        matrix = ref_cell.dot(np.linalg.inv(strained_cell))
        identity = np.identity(3)
        strain = matrix - identity

        # Verify that the strain matrix is symmetric
        assert np.allclose(strain, strain.T)

        if not principal:
            return strain

        return np.linalg.eigvalsh(strain)
