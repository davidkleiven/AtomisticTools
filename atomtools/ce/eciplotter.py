import numpy as np
from matplotlib import pyplot as plt

class ECIPlotter(object):
    def __init__(self, eci, naming="raw"):
        allowed_namings = ["raw", "normalized"]
        self.eci = eci
        if naming not in allowed_namings:
            raise ValueError("Naming has to be one of {}".format(
                allowed_namings))
        self.naming = naming
        self._rotation = {
            "raw": "vertical",
            "normalized": "horizontal"
        }

    def sort_on_cluster_size( self ):
        """
        Sort the atoms according to cluster size
        """
        cluster_size = []
        extents = {}
        eci_names = {}
        eci_values = {}

        for key,value in self.eci.items():
            size = int(key[1])
            extent =  cluster_dia_from_name(key)
            if ( size in eci_names.keys() ):
                eci_names[size].append(key)
                eci_values[size].append(value)
                extents[size].append(extent)
            else:
                eci_names[size] = [key]
                eci_values[size] = [value]
                extents[size] = [extent]

        # Within each size sort on absolute value
        for key,value in eci_values.items():
            indx_srt = np.argsort(np.abs(extents[key]))
            new_list = [value[indx] for indx in indx_srt]
            eci_values[key] = new_list
            new_list = [eci_names[key][indx] for indx in indx_srt]
            eci_names[key] = new_list
        return eci_names, eci_values

    def plot( self, show_names=True ):

        labels = {
            0:"Bias",
            1:"Singlets",
            2:"Doublets",
            3:"Triplets",
            4:"Quadruplets"
        }
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        eci_names, eci_vals = self.sort_on_cluster_size()
        prev = 0
        all_keys = []
        all_indx = []
        for key, value in eci_vals.items():
            all_keys += eci_names[key]
            indx = np.arange(prev,prev+len(value))
            all_indx += list(indx)
            ax.bar(indx, value, label=labels[key])
            prev = indx[-1]+2
        ax.legend( loc="best", frameon=False )
        ax.axhline(0.0, color="black", linewidth=0.5, ls="--")
        ax.set_xticklabels([])
        ax.set_ylabel( "ECI (eV/atom)" )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if show_names:
            ax.set_xticks(all_indx)
            names = self._get_names(all_keys)
            ax.set_xticklabels(names, rotation=self._rotation[self.naming])
            if self.naming == "normalized":
                ax.set_xlabel("$d/r_{nn}")
        return fig

    def _get_names(self, raw_names):
        """Return proper names."""
        if self.naming == "raw":
            return raw_names
        elif self.naming == "normalized":
            diameters = []
            names = []
            for name in raw_names:
                dia = cluster_dia_from_name(name)
                if dia > 0.01:
                    diameters.append(dia)
                else:
                    names.append(0)
            min_dia = np.min(diameters)
            names += [np.round(dia/min_dia, 2) for dia in diameters]
            return names


def cluster_dia_from_name( cname ):
    if ( int(cname[1]) <= 1 ):
        return 0.0

    splitted = cname.split("_")
    if ( splitted[1].find("p") != -1 ):
        float_version = float( splitted[1].replace("p",".") )
        return int(float_version*100)

    return int(splitted[1])
