from ase.db import connect
import numpy as np

class BondLengthDistribution(object):
    """Class for exploring how much bond lengths change."""
    def __init__(self, atoms, num_bonds=2):
        self.num_bonds = num_bonds
        self.atoms = atoms
        self.bond_lengths = {}

    @staticmethod
    def symbs2key(symbs):
        """Convert tuple of two symbols into a key."""
        symbs = sorted(symbs)
        return "-".join(symbs)

    def _get_closest(self, ref_indx):
        """Find the atoms being closest."""
        indices = list(range(0, len(self.atoms)))
        indices.remove(ref_indx)
        dists = self.atoms.get_distances(ref_indx, indices, mic=True)
        srt_indx = np.argsort(dists)
        closest_indices = [indices[srt_indx[i]] for i in range(self.num_bonds)]

        # Remove indices smaller than ref_indx
        filtered_closest_indx = []
        for indx in closest_indices:
            if indx > ref_indx:
                filtered_closest_indx.append(indx)
        closest_indices = filtered_closest_indx
        closest_dists = [dists[srt_indx[i]] for i in range(self.num_bonds)]
        symbs = [(self.atoms[ref_indx].symbol, self.atoms[indx].symbol) for indx in closest_indices]
        
        # Create a dictionary with the symbols and distances
        dist_dict = {}
        for dist, symb in zip(closest_dists, symbs):
            key = BondLengthDistribution.symbs2key(symb)
            if key not in dist_dict.keys():
                dist_dict[key] = [dist]
            else:
                dist_dict[key].append(dist)
        return dist_dict

    def calculate_bond_lengths(self):
        """Calculate bond length changes."""
        self.bond_lengths = {}
        for ref_indx in range(len(self.atoms)):
            dist_dict = self._get_closest(ref_indx)
            for k, v in dist_dict.items():
                if k not in self.bond_lengths.keys():
                    self.bond_lengths[k] = v
                else:
                    self.bond_lengths[k] += v
        return self.bond_lengths

    def show_bond_lengths(self, show=False):
        """Create plots of the bond lengths."""
        from matplotlib import pyplot as plt
        num_figs = len(self.bond_lengths.keys())
        if num_figs < 3:
            num_cols = num_figs
            num_rows = 1
        else:
            num_cols = 3
            num_rows = int(num_figs/num_cols)
            if num_rows < float(num_figs)/num_cols:
                num_rows += 1
        fig, ax = plt.subplots(num_rows, num_cols)

        ax_count = 0
        for k, v in self.bond_lengths.items():
            if isinstance(ax, np.ndarray):
                col = ax_count%num_cols
                row = int(ax_count/num_cols)
                current_ax = ax[row, col]
            else:
                current_ax = ax
            current_ax.plot(v, "x")
            current_ax.set_title(k)
            mean = np.mean(v)
            current_ax.axhline(mean, ls="--")
            ax_count += 1
        
        if show:
            plt.show()
        return fig


