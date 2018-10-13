import unittest
import os

db_name = "test_db.db"
class TestPhononDOS( unittest.TestCase ):
    def test_phonon_dos(self):
        no_throw = True
        print ( "NOTE! This test does not make sense at the moment as CE is not part of the ASE distribution available via pip")
        try:
            from atomtools.ce import phonon_ce_eval as pce
            from ase.clease.settings import CEBulk
            conc_args = {
                "conc_ratio_min_1":[[1,0]],
                "conc_ratio_max_1":[[0,1]],
            }
            BC = CEBulk( "fcc", 4.05, None, [4,4,4], 1, [["Al","Mg"]], conc_args, db_name, max_cluster_size=4, reconf_db=False)
            ph = pce.PhononEval( BC )
        except ImportError as exc:
            print ("Import error! This expected due to CE missing on pip")
        except Exception as exc:
            print ( str(exc) )
            no_throw = False
        self.assertTrue( no_throw )

    def __del__(self):
        if ( os.path.isfile(db_name) ):
            os.remove(db_name)

if __name__ == "__main__":
    unittest.main()
