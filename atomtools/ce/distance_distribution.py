import numpy as np

class DistanceDistribution(object):
    """Calculate the distribution of distances in an object."""
    def __init__(self, atoms, exclude=["X"]):
        self.atoms = atoms
        self.exclude = exclude
        self._distances = self._calculate_distances()
    
    def _calculate_distances(self):
        """Calculates all internal distances."""
        all_dists = []
        for ref in range(len(self.atoms)):
            if self.atoms[ref].symbol in self.exclude:
                continue
            indices = list(range(ref+1, len(self.atoms)))
            indices = self._filter_excluded(indices)
            if len(indices) == 0:
                continue
            dists = self.atoms.get_distances(ref, indices, mic=True)
            all_dists += list(dists)
        
        # Normalize by the mean distance
        return np.array(all_dists)/np.mean(all_dists)

    def _filter_excluded(self, indices):
        """Remove excluded indices."""
        return [indx for indx in indices 
                if self.atoms[indx].symbol not in self.exclude]

    def get_cumulative_distribution(self):
        """Return the cumulative distribution."""
        srt_dists = np.sort(self._distances)
        tol = 1E-3
        for i in range(1, len(srt_dists)):
            while srt_dists[i] - srt_dists[i-1] < tol:
                srt_dists[i] += tol
        return {"x": srt_dists-srt_dists[0], 
                "P": np.linspace(0.0, 1.0, len(self._distances), endpoint=False)}

    def kolmogorov_smirnov_distance(self, other, show=False):
        """Calculates the Kolmogorov-Smirnov distance.
        
        :param DistanceDistribution other: Other distance distribution method
        """
        from scipy.interpolate import interp1d
        cumdist1 = self.get_cumulative_distribution()
        cumdist2 = other.get_cumulative_distribution()
       
        # Normalize the x-values
        # range1 = np.max(cumdist1["x"]) - np.min(cumdist1["x"])
        # range2 = np.max(cumdist2["x"]) - np.min(cumdist2["x"])
        # cumdist1["x"] -= np.min(cumdist1["x"])
        # cumdist1["x"] *= (range2/range1) 
        # cumdist1["x"] += np.min(cumdist2["x"])

        interp_cumdist1 = interp1d(cumdist1["x"], cumdist1["P"], kind="linear",
                                   fill_value=(0.0, 1.0), bounds_error=False)

        diff = cumdist2["P"] - interp_cumdist1(cumdist2["x"])

        ks_distance = np.max(np.abs(diff))

        if show:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(cumdist1["x"], cumdist1["P"], drawstyle="steps", label="First")
            ax.plot(cumdist2["x"], cumdist2["P"], drawstyle="steps", label="Second")
            ax.set_xlabel("Normalized distances")
            ax.set_ylabel("Cummulative Distribution")
            ax.legend(loc="best")
            plt.show()
        return ks_distance


            