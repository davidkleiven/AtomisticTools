import numpy as np


class GaussianMixtureClassifier(object):
    """Class for classifying data points via Gaussian Mixture.
    :param evaluator: Instance of the Evaluator class in ASE
    """

    def __init__(self, evaluator):
        from ase.ce import Evaluate
        from copy import deepcopy
        if not isinstance(evaluator, Evaluate):
            raise TypeError("evaluator has to be of type Evaluate")

        self.evaluator = evaluator
        self.orig_cf_matrox = deepcopy(self.evaluator.cf_matrix)

    def classify(self, alpha=1E-5):
        from sklearn.mixture import GaussianMixture
        E_dft = self.evaluator.e_dft
        self.evaluator.get_eci(alpha)
        E_pred = self.evaluator.cf_matrix.dot(self.evaluator.eci)
        self.evaluator.cv_loo(alpha)
        mix_model = GaussianMixture(n_components=3, init_params='random',
                                    n_init=10)
        data = np.column_stack((E_pred, E_dft))
        mix_model.fit(data)
        labels = mix_model.predict(data)
        self.show(data, labels)

    def show(self, data, labels):
        """Visualize the results of the mixture model."""
        from matplotlib import pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        max_label = np.max(labels)
        for group in range(max_label+1):
            x = data[:, 0][labels == group]
            y = data[:, 1][labels == group]
            ax.plot(x, y, "o", mfc="none")
        return fig
