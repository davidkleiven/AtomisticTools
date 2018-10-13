import unittest
import os

db_name = "test_db.db"
class TestPopulationCovariance( unittest.TestCase ):
    def test_pop_cov(self):
        no_throw = True
        try:
            from atomtools.ce import population_variance as pv
            from ase.clease.settings import CEBulk
            conc_args = {
                "conc_ratio_min_1":[[1,0]],
                "conc_ratio_max_1":[[0,1]],
            }
            BC = CEBulk( "fcc", 4.05, None, [4,4,4], 1, [["Al","Mg"]], conc_args, db_name, max_cluster_size=4, reconf_db=False)
            pvcov = pv.PopulationVariance(BC)
        except ImportError as exc:
            print (str(exc))
            print ("ImportError! Expected as ASE on pip has no CE")
        except Exception as exc:
            no_throw = False
            print (str(exc))
        self.assertTrue(no_throw)

    def __del__(self):
        if ( os.path.isfile(db_name) ):
            os.remove(db_name)

if __name__ == "__main__":
    unittest.main()
