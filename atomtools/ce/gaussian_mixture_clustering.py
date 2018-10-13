import numpy as np
import dataset


class GaussianMixtureClassifier(object):
    """Class for classifying data points via Gaussian Mixture.
    :param evaluator: Instance of the Evaluator class in ASE
    """

    def __init__(self, evaluator):
        from ase.clease import Evaluate
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


class FilterCollapsed(object):
    """Class that attempts to filter out collapsed structures."""

    def __init__(self, evaluator, db_name="filter_db.db", restart=True):
        from ase.clease import Evaluate
        from copy import deepcopy
        if not isinstance(evaluator, Evaluate):
            raise TypeError("evaluator has to be of type Evaluate")
        self.evaluator = evaluator
        self.orig_cf_matrix = deepcopy(evaluator.cf_matrix)
        self.orig_e_dft = deepcopy(self.evaluator.e_dft)
        self.names_removed = []
        self.db_name = db_name
        if restart:
            self.remove_already_calculated()

    def filter_worst(self, alpha):
        """Remove the point that is furthest away."""
        E_dft = self.evaluator.e_dft
        self.evaluator.get_eci(alpha)
        E_pred = self.evaluator.cf_matrix.dot(self.evaluator.eci)
        diff = E_dft - E_pred
        min_indx = np.argmin(diff)
        name_removed = self._remove_indx(min_indx)
        return name_removed

    def _remove_indx(self, indx):
        name_removed = self.evaluator.names[indx]
        self.evaluator.e_dft = np.delete(self.evaluator.e_dft, indx)
        del self.evaluator.names[indx]
        self.evaluator.cf_matrix = \
            np.delete(self.evaluator.cf_matrix, indx, axis=0)
        return name_removed

    def remove_already_calculated(self):
        """Remove already calculated."""
        db = dataset.connect("sqlite:///{}".format(self.db_name))
        tbl = db["unique_names"]
        for entry in tbl.find():
            name = entry["name"]
            indx = self.evaluator.names.index(name)
            self._remove_indx(indx)
            self.names_removed.append(name)

    def _find_new_best_alpha(self, alpha_min, alpha_max, num_alpha):
        """Find the alpha value that is best."""
        return self.evaluator.plot_CV(alpha_min, alpha_max,
                                      num_alpha=num_alpha)

    def run(self, alpha, n_points=10, update_alpha=True,
            alpha_min=1E-5, alpha_max=1E-2, num_alpha=8):
        """Filter until converged."""
        db = dataset.connect("sqlite:///{}".format(self.db_name))
        tbl_stat = db["status"]
        tbl_unique_name = db["unique_names"]
        for _ in range(n_points):
            name = self.filter_worst(alpha)
            cv = self.evaluator.cv_loo(alpha) * 1000.0
            row = {"alpha": alpha, "cv": cv, "name": name}
            tbl_stat.insert(row)
            if name not in self.names_removed:
                tbl_unique_name.insert({"name": name})
            print("Removed: {}. New CV: {} meV/atom".format(name, cv))
            if update_alpha:
                alpha = self._find_new_best_alpha(alpha_min, alpha_max,
                                                  num_alpha)

    @staticmethod
    def get_selection_condition(db_name, cv_score):
        """Return a selection criteria based on the values to be left out."""
        db = dataset.connect("sqlite:///{}".format(db_name))
        tbl_name = db["unique_names"]
        tbl_stat = db["status"]
        scond = []
        for row in tbl_name.find():
            cv = tbl_stat.find_one(name=row["name"])["cv"]
            if cv > cv_score:
                scond.append(("name", "!=", row["name"]))
        return scond
