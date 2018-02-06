import unittest

class TestPopulationCovariance( unittest.TestCase ):
    def test_pop_cov(self):
        no_throw = True
        try:
            from atomtools.ce import population_variance as pv
            BC = None
            pvcov = pv.PopulationVariance(BC)
        except ImportError as exc:
            print (str(exc))
            print ("ImportError! Expected as ASE on pip has no CE")
        except Exception as exc:
            no_throw = False
            print (str(exc))
        self.assertTrue(no_throw)

if __name__ == "__main__":
    unittest.main()
