import numpy as np


class AutocorrelationFunction(object):
    """
    Calculate the auto correlation function of an atoms object
    along certain direction.
    """

    def __init__(self, atoms=None, plane_normal=(1, 1, 1), symb_dict={}):
        self.atoms = atoms
        self.plane_normal = np.array(plane_normal)
        self.symb_dict = symb_dict

    @property
    def projected_positions(self):
        """Project the positions onto a plane."""
        normal_cart = self.plane_normal / \
            np.sqrt(self.plane_normal.dot(self.plane_normal))
        pos = self.atoms.get_positions()

        # Rotate such that z axis points along the normal vector
        alpha = np.arctan2(normal_cart[1], normal_cart[0])
        sa = np.sin(alpha)
        ca = np.cos(alpha)
        R = np.eye(3)
        R[0, 0] = ca
        R[0, 1] = -sa
        R[1, 0] = sa
        R[1, 1] = ca
        pos = (R.dot(pos.T)).T
        print(R)

        alpha = np.arccos(normal_cart[2])
        sa = np.sin(alpha)
        ca = np.cos(alpha)
        R = np.eye(3)
        R[0, 0] = ca
        R[0, 2] = -sa
        R[2, 0] = sa
        R[2, 2] = ca
        pos = (R.dot(pos.T)).T
        # assert np.allclose(pos[:, 2], np.zeros(pos.shape[0]))
        return pos[:, :2]

    @property
    def unique_symbs(self):
        un_symb = []
        for atom in self.atoms:
            if atom.symbol not in un_symb:
                un_symb.append(atom.symbol)
        return un_symb

    def projected_image(self, npix=512, show=True, cmap="gray"):
        """Construct a projected image."""
        from scipy.signal import correlate2d
        for symb in self.unique_symbs:
            if symb not in self.symb_dict.keys():
                self.symb_dict[symb] = 0

        values = np.zeros(len(self.atoms))
        for i, atom in enumerate(self.atoms):
            values[i] = self.symb_dict[atom.symbol]
        pos = self.projected_positions

        image = np.zeros((npix, npix))
        xmin = np.min(pos[:, 0])
        xmax = np.max(pos[:, 0])
        ymin = np.min(pos[:, 1])
        ymax = np.max(pos[:, 1])
        dx = (xmax - xmin)/(npix-1)
        dy = (ymax - ymin)/(npix-1)

        for i in range(pos.shape[0]):
            ix = int((pos[i, 0] - xmin)/dx)
            iy = int((pos[i, 1] - ymin)/dy)
            image[ix, iy] += values[i]
        corr = correlate2d(image, image)

        if show:
            self.plot(image, cmap=cmap)
            self.plot(corr, cmap=cmap)
        return image, corr

    def plot(self, image, cmap="grey"):
        """Create a plot of the image and the show the correlation function."""
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(image, cmap=cmap)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        return fig
