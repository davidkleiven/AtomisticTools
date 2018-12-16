import numpy as np

class FilterDisplacements(object):
    def __init__(self, max_displacement=0.2, exclude=["X"]):
        self.max_displacement = max_displacement
        self.exclude = exclude
        self.rejected_reason = ""

    def is_valid(self, atoms_init, atoms_final, num_attempts=1):
        for _ in range(num_attempts):
            if self.is_valid_single_attempt(atoms_init, atoms_final):
                return True
        return False

    def is_valid_single_attempt(self, atoms_init, atoms_final):
        """Valid if the maximum displacement is less than threshold."""
        from scipy.spatial import cKDTree as KDTree
        from random import shuffle
        atoms1 = atoms_init.copy()
        atoms2 = atoms_final.copy()

        vol1 = atoms1.get_volume()
        vol2 = atoms2.get_volume()
        if vol2 > vol1:
            ratio = (vol2/vol1)**(1.0/3.0)
            cell1 = atoms1.get_cell()
            atoms1.set_cell(cell1*ratio, scale_atoms=True)
        else:
            ratio = (vol1/vol2)**(1.0/3.0)
            cell2 = atoms2.get_cell()
            atoms2.set_cell(cell2*ratio, scale_atoms=True)

        # Try construct the relation
        used_indices = []
        tree = KDTree(atoms2.get_positions())
        indices = list(range(0, len(atoms1)))
        shuffle(indices)
        for atom in atoms1:
            if atom.symbol in self.exclude:
                continue
            dist, closest = tree.query(atom.position, k=12)
            srt_indx = np.argsort(dist)
            dist = [dist[indx] for indx in srt_indx]
            closest = [closest[indx] for indx in srt_indx]

            if all(c in used_indices for c in closest):
                # More than one atom is closest to this
                # structure
                self.rejected_reason = "More than one atom mapped onto the "
                self.rejected_reason += "same atoms in the initial structure"
                return False

            # First, unused with mathing symbol
            closest_indx = None
            closest_dist = None
            for i, indx in enumerate(closest):
                if atoms2[indx].symbol == atom.symbol and indx not in used_indices:
                    closest_indx = indx
                    closest_dist = dist[i]
                    break

            if closest_indx is None:
                self.rejected_reason = "No unused atoms with macthing symbol!"
                return False
            
            used_indices.append(closest_indx)
            if closest_dist > self.max_displacement:
                # The displacement is larger than the tolereance
                self.rejected_reason = "Max displacement too large"
                return False
            
            if atom.symbol != atoms2[closest_indx].symbol:
                self.rejected_reason = "Mapped symbol does not match!"
                return False
        return True
            
            