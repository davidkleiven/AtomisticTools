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

    def _to_mandel(self, tensor):
        """Convert tensor to mandel notation."""
        mandel = np.zeros(6)
        mandel[0] = tensor[0, 0]
        mandel[1] = tensor[1, 1]
        mandel[2] = tensor[2, 2]
        mandel[3] = np.sqrt(2.0)*tensor[1, 2]
        mandel[4] = np.sqrt(2.0)*tensor[0, 2]
        mandel[5] = np.sqrt(2.0)*tensor[0, 1]
        return mandel

    def _to_mandel_rank4(self, tensor):
        """Convert rank 4 tensor to mandel notation."""
        from itertools import product
        out = np.zeros((6, 6))
        mandel_lut = {
            (0, 0): 0,
            (1, 1): 1,
            (2, 2): 2,
            (1, 2): 3,
            (0, 2): 4,
            (0, 1): 5
        }
        for ind in product([0, 1, 2], repeat=4):
            if ind[1] < ind[0] or ind[3] < ind[2]:
                continue
            row = mandel_lut[(ind[0], ind[1])]
            col = mandel_lut[(ind[2], ind[3])]

            if ind[0] != ind[1] and ind[2] != ind[3]:
                value = 2.0*tensor[ind[0], ind[1], ind[2], ind[3]]
            elif ind[0] != ind[1]:
                value = np.sqrt(2.0)*tensor[ind[0], ind[1], ind[2], ind[3]]
            elif ind[2] != ind[3]:
                value = np.sqrt(2.0)*tensor[ind[0], ind[1], ind[2], ind[3]]
            else:
                value = tensor[ind[0], ind[1], ind[2], ind[3]]
            out[row, col] = value
        return out

    def _to_full_rank4(self, mandel_tensor):
        """Convert Mandel representation to full tensor."""
        from itertools import product
        out = np.zeros((3, 3, 3, 3))
        mandel_lut = [(0, 0), (1, 1), (2, 2), 
                      (1, 2), (0, 2), (0, 1)]
        
        for ind in product(range(6), repeat=2):
            if ind[0] > 2 and ind[1] > 2:
                value = mandel_tensor[ind[0], ind[1]]/2.0
            elif ind[0] > 2 or ind[1] > 2:
                value = mandel_tensor[ind[0], ind[1]]/np.sqrt(2.0)
            else:
                value = mandel_tensor[ind[0], ind[1]]

            row = mandel_lut[ind[0]]
            col = mandel_lut[ind[1]]

            out[row[0], row[1], col[0], col[1]] = value
            out[row[0], row[1], col[1], col[0]] = value
            out[row[1], row[0], col[1], col[0]] = value
            out[row[1], row[0], col[0], col[1]] = value
        return out

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
                strain = self._to_mandel(strain)
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
                strain = self._to_mandel(strain)
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

    def get(self, select_cond=[], strains=None, stresses=None, spg=1, perm="xyz",
            convert_stress_to_mandel=True, convert_strain_to_mandel=False):
        """Compute the elastic properties."""
        if strains is not None and stresses is not None:
            for eps, sigma in zip(strains, stresses):
                d = {
                    "strain": eps,
                    "stress": sigma
                }
                self.data.append(d)
        else:
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

        # Convert to Mandel stress (strain is already assumed to be
        # in mandel representation)
        if convert_stress_to_mandel:
            stress_matrix[3:, :] *= np.sqrt(2.0)

        if convert_strain_to_mandel:
            strain_matrix[3:, :] *= np.sqrt(2.0)
        print("here")
        print(stress_matrix)
        print("here")
        prec = np.linalg.inv(strain_matrix.dot(strain_matrix.T))
        #self.elastic_tensor = prec.dot(stress_matrix.dot(strain_matrix.T))
        self.elastic_tensor = stress_matrix.dot(strain_matrix.T).dot(prec)
        self._symmetrize_elastic_tensor(spg=spg, perm=perm)
        return self.elastic_tensor

    def _symmetrize_elastic_tensor(self, spg=1, perm="xyz"):
        if spg == 1:
            return
        permut_lut = {
            "xyz": 0,
            "zxy": -1,
            "yzx": -2
        }
        from ase.spacegroup import Spacegroup
        spg = Spacegroup(spg)
        sym_op = spg.get_rotations()
        full = self._to_full_rank4(self.elastic_tensor)
        new_tensor = np.zeros((3, 3, 3, 3))
        for op in sym_op:
            op = np.roll(op, permut_lut[perm], (0, 1))
            avg_tensor = np.zeros((3, 3, 3, 3))
            avg_tensor = np.einsum("pl,ijkl->ijkp", op, full)
            avg_tensor = np.einsum("ok,ijkp->ijop", op, avg_tensor)
            avg_tensor = np.einsum("nj,ijop->inop", op, avg_tensor)
            avg_tensor = np.einsum("mi,inop->mnop", op, avg_tensor)
            new_tensor += avg_tensor
        new_tensor /= len(sym_op)
        self.elastic_tensor = self._to_mandel_rank4(new_tensor)

    @property
    def compliance_tensor(self):
        """Return the compliance tensor."""
        if self.elastic_tensor is None:
            msg = "Elastic tensor not computed."
            msg += "Call get() method first."
            raise ValueError(msg)
        
        # Convert Mandel tensor into Voigt tensor
        C = self.elastic_tensor.copy()
        C[3:, :3] /= np.sqrt(2.0)
        C[:3, 3:] /= np.sqrt(2.0)
        C[3:, 3:] /= 2.0
        return np.linalg.inv(C)

    @property
    def isotropic_elastic_tensor(self):
        """Return the isotropic elastic tensor.
        Uses the bulk modulus and the poisson ratio obtained from the
        general tensor (in Mandel notation)
        """
        K = self._bulk_mod_voigt()
        mu = self._shear_mod_voigt()
        tensor = np.zeros((6, 6))
        tensor[0, 0] = tensor[1, 1] = tensor[2, 2] = K + 4.0*mu/3.0
        tensor[3, 3] = tensor[4, 4] = tensor[5, 5] = 2*mu
        tensor[0, 1] = tensor[0, 2] = K - 2.0 * mu/3.0
        tensor[1, 0] = tensor[1, 2] = K - 2.0 * mu/3.0
        tensor[2, 0] = tensor[2, 1] = K - 2.0 * mu/3.0
        return tensor

    def _bulk_mod_voigt(self):
        """Bulk modulus by the voigt average."""
        if self.elastic_tensor is None:
            msg = "Elastic tensor not computed."
            msg += "Call get() method first."
            raise ValueError(msg)
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
        """Shear modulus by voigt average."""
        if self.elastic_tensor is None:
            msg = "Elastic tensor not computed."
            msg += "Call get() method first."
            raise ValueError(msg)
        C = self.elastic_tensor
        Gv = C[0, 0] + C[1, 1] + C[2, 2]
        Gv -= (C[0, 1] + C[1, 2] + C[2, 0])
        Gv += 3.0*(C[3, 3] + C[4, 4] + C[5, 5])*0.5
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
