import unittest
from atomtools.eos.birch_murnagan import BirschMurnagan
import numpy as np
V = np.linspace( 10.0,20.0, 100 )
E = 2.0 + np.exp(-V/15.0) # Some random function decaying to a constant
class TestEOS( unittest.TestCase ):
    def test_bm_eos(self):
        no_throw = True
        msg = ""
        try:
            bm = BirschMurnagan( V, E )
            bm.fit()
            vol = np.linspace(5.0,30.0,50 )
            B = bm.bulk_modulus( vol )
        except Exception as exc:
            msg = str(exc)
            no_throw = False
        self.assertTrue( no_throw, msg=msg )

if __name__ == "__main__":
    unittest.main()
