import unittest
import matplotlib as mpl
mpl.use("Agg") # This does not require any screen
msg = ""
try:
    from atomtools.ce.eciplotter import ECIPlotter
    available = True
except ImportError as exc:
    available = False
    msg = str(exc)

ecis = {
    "c1_1":1.0,
    "c2_1000_1":-1.0,
    "c2_1225_1":3.0,
    "c3_707_1":-4.0
}
class TestECIPlotter(unittest.TestCase):
    def test_plot(self):
        if ( not available ):
            self.skipTest( "Test not available: {}".format(msg) )
            return
        no_throw = True
        try:
            plotter = ECIPlotter(ecis)
            plotter.plot()
        except Exception as exc:
            print (str(exc))
            no_throw = False
        self.assertTrue( no_throw )

if __name__ == "__main__":
    unittest.main()
