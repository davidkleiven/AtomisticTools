from ase.db import connect
import numpy as np

def plot_normalized_bond_lengths(atoms_list, num_bonds=2, plot_args=[]):
    """Create a plot with bond lengths of list.
    
    :param atoms_list: Nested list of atoms object
    :type atoms_list: list of lists of Atoms objects
    :param int num_bonds: Number of nearest neighbours to consider
    :param plot_args: List of dictionary with plotting arguments
        has to be the same length as atom_list
        Example: {"alpha": 0.5, "marker": "s"}
    """
    from matplotlib import pyplot as plt
    data = []
    for i, atom_grp in enumerate(atoms_list):
        data.append({})
        for atom in atom_grp:
            b_length = BondLengthDistribution(atom, num_bonds=num_bonds)
            new_data = b_length.calculate_bond_lengths(normalize_by_mean=True)
            for k, v in new_data:
                if k not in data[i].keys():
                    data[k] = v
                else:
                    data[k] += v

    assert len(plot_args) == len(atoms_list)
    all_keys = []
    for dset in data:
        all_keys.append(list(dset.keys()))
    all_keys = sorted(list(set(all_keys)))

    num_figs = len(all_keys)
    if num_figs < 3:
        num_cols = num_figs
        num_rows = 1
    else:
        num_cols = 3
        num_rows = int(num_figs/num_cols)
        if num_rows < float(num_figs)/num_cols:
            num_rows += 1
    fig, ax = plt.subplots(num_rows, num_cols)
    
    for plot_args, dset in zip(data, plot_args):
        ax_count = 0
        for k in all_keys:
            if k not in dset.keys():
                continue
        if isinstance(ax, np.ndarray):
            col = ax_count%num_cols
            row = int(ax_count/num_cols)
            current_ax = ax[row, col]
        else:
            current_ax = ax
        current_ax.plot(dset[k], **plot_args)
    return fig


def bond_lengths_score(init_struct, final_struct, num_bonds=2, show=False):
    """Calculates a score based on the bond lengths."""    
    b_init = BondLengthDistribution(init_struct, num_bonds=num_bonds)
    b_lengths_init = b_init.calculate_bond_lengths(normalize_by_mean=True)

    b_final = BondLengthDistribution(final_struct, num_bonds=num_bonds)
    b_lengths_final = b_final.calculate_bond_lengths(normalize_by_mean=True)

    # Calculate the cummulative distribution of the bond lengths
    cumdist_init = {}
    for k, v in b_lengths_init.items():
        if np.allclose(v, 1.0):
            raise ValueError("There is only one bond type in the initial structure. "
                             "Maybe increase the number of bonds.")
        cumdist_init[k] = {"x": np.sort(v), 
                           "P": np.linspace(0.0, 1.0, len(v), endpoint=False)}
    cumdist_final = {}
    for k, v in b_lengths_final.items():
        cumdist_final[k] = {"x": np.sort(v), 
                            "P": np.linspace(0.0, 1.0, len(v), endpoint=False)}

    colors = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9']
    keys = list(cumdist_init.keys()) + list(cumdist_final.keys())
    keys = list(set(keys))
    if show:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        count = 0
        for k in keys:
            if k in cumdist_init.keys() and k in cumdist_final.keys():
                ax.plot(cumdist_init[k]["x"], cumdist_init[k]["P"], 
                        linestyle="-", drawstyle="steps",
                        color=colors[count%len(colors)])

                ax.plot(cumdist_final[k]["x"], cumdist_final[k]["P"], 
                        linestyle="--", drawstyle="steps", 
                        color=colors[count%len(colors)])
                count += 1
        plt.show()


class BondLengthDistribution(object):
    """Class for exploring how much bond lengths change."""
    def __init__(self, atoms, num_bonds=2, exclude=["X"]):
        self.num_bonds = num_bonds
        self.atoms = atoms.copy()
        self.bond_lengths = {}
        self.exclude_symbs = exclude

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
        closest_indices = []
        
        num_inserted = 0
        for indx in srt_indx:
            if self.atoms[indx].symbol in self.exclude_symbs:
                continue
            closest_indices.append(indx)
            num_inserted += 1
            if num_inserted >= self.num_bonds:
                break
        
        if len(closest_indices) != self.num_bonds:
            raise RuntimeError("Could not find the requested number of neigbours!")

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

    def calculate_bond_lengths(self, normalize_by_mean=False):
        """Calculate bond length changes."""
        self.bond_lengths = {}
        for ref_indx in range(len(self.atoms)):
            dist_dict = self._get_closest(ref_indx)
            for k, v in dist_dict.items():
                if k not in self.bond_lengths.keys():
                    self.bond_lengths[k] = v
                else:
                    self.bond_lengths[k] += v
        
        if normalize_by_mean:
            for k, v in self.bond_lengths.items():
                self.bond_lengths[k] /= np.mean(v)
        return self.bond_lengths

    def show_bond_lengths(self, show=False, fig=None):
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


