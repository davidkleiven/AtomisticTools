import unittest

class TestPhononDOS( unittest.TestCase ):
    def test_phonon_dos(self):
        no_throw = True
        print ( "NOTE! This test does not make sense at the moment as CE is not part of the ASE distribution available via pip")
        try:
            from atomtools.ce import phonon_ce_eval as pce
            BC = None # This should be a bulk crystal object
            ph = pce.PhononEval( BC )
        except ImportError as exc:
            print ("Import error! This expected due to CE missing on pip")
        except Exception as exc:
            print ( str(exc) )
            no_throw = False
        self.assertTrue( no_throw )

if __name__ == "__main__":
    unittest.main()
