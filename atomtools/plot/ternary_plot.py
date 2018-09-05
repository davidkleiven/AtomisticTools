from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


class TernaryPlot(object):
    def __init__(self, x, y, z, labels=["x", "y", "z"], color=None,
                 cbar_label="Colorbar"):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        self.color = color
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.labels = labels
        self.cbar_label = cbar_label

    def _create_vertices(self, lw=2, num_ticks=4, grid_color="#bdbdbd",
                         grid_lw=0.5):
        """Create the vertices of the triangle."""
        # Draw grid
        tick_points = np.linspace(0.0, 1.0, 5)
        S = self.x[0] + self.y[0] + self.z[0]
        for x0 in tick_points:
            x_cart = [0.5 * (S - x0) / S, (S - x0) / S]
            y_cart = [0.5 * np.sqrt(3) * (S-x0)/S, 0]
            self.ax.plot(x_cart, y_cart, color=grid_color, lw=grid_lw)
            self.ax.annotate(xy=(x_cart[1], y_cart[1]), s=str(x0),
                             xytext=(x_cart[1], y_cart[1] - 0.07))

        for y0 in tick_points:
            x_cart = [y0 / S, 0.5 * (S + y0)/S]
            y_cart = [0, 0.5 * np.sqrt(3) * (S - y0) / S]
            self.ax.plot(x_cart, y_cart, color=grid_color, lw=grid_lw)
            self.ax.annotate(xy=(x_cart[1], y_cart[1]), s=str(y0),
                             xytext=(x_cart[1] + 0.0, y_cart[1]))

        for z0 in tick_points:
            x_cart = [0.5 * z0 / S, 0.5 * (2 * S - z0) / S]
            y_cart = [0.5 * np.sqrt(3) * z0 / S, 0.5 * np.sqrt(3) * z0 / S]
            self.ax.plot(x_cart, y_cart, color=grid_color, lw=grid_lw)
            self.ax.annotate(xy=(x_cart[1], y_cart[1]), s=str(z0),
                             xytext=(x_cart[0] - 0.125, y_cart[0]))

        self.ax.plot([0, 1], [0, 0], color="black")
        self.ax.plot([0, 0.5], [0, np.sqrt(3)/2.0], color="black")
        self.ax.plot([1, 0.5], [0, np.sqrt(3)/2.0], color="black")
        self.ax.text(0.5, -0.1, self.labels[0])
        self.ax.text(0.75, 0.5, self.labels[1])
        self.ax.text(0.2, 0.5, self.labels[2])


    def plot(self, cmap="copper", **plt_args):
        """Create a ternay plot."""
        self._create_vertices()

        S = self.x + self.y + self.z
        x_cartesian = 0.5 * (2*self.y + self.z) / S
        y_cartesian = 0.5 * np.sqrt(3) * self.z / S

        if "marker" not in plt_args.keys():
            plt_args["marker"] = "v"

        if self.color is not None:
            norm = Normalize(vmin=np.min(self.color), vmax=np.max(self.color))
            smap = ScalarMappable(norm=norm, cmap=cmap)
            color = [smap.to_rgba(c) for c in self.color]
            plt_args["c"] = color
        im = self.ax.scatter(x_cartesian, y_cartesian, **plt_args)

        # Remove spines
        self.ax.spines["right"].set_visible(False)
        self.ax.spines["top"].set_visible(False)
        self.ax.spines["left"].set_visible(False)
        self.ax.spines["bottom"].set_visible(False)
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticklabels([])
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        if self.color is not None:
            ticks = np.linspace(np.min(self.color), np.max(self.color), 4)
            smap.set_array([])
            cbar = self.fig.colorbar(smap, ax=self.ax, label=self.cbar_label)
