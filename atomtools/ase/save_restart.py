
class SaveRestartFiles(object):
    def __init__(self, calc, name):
        """
        Callable class that can be passed attached to ASE relaxers

        Arguments
        ----------
        calc - Instance of the calculator
        name - Unique name that will be included in the filename
        """
        self.calc = calc
        self.name = name

    @staticmethod
    def restart_name(name):
        """
        Returns the filename that will be used in the simulations
        """
        return fname = "calc_restart{}.gpw".format(name)

    def __call__(self):
        """
        Stores the current state of the calculator by invoking its write
        method
        """
        fname = SaveRestartFiles.restart_name(self.name)
        self.calc.write(fname)
