import unittest
msg = ""
try:
    from atomtools.ce import ce_phonon_dos as cpd
    available = True
except ImportError as exc:
    available = False
    msg = str(exc)
import numpy as np
import os

db_name = "test_phonon.db"
class TestPhononDB(unittest.TestCase):
    def test_save_load(self):
        if ( not available ):
            self.skipTest("Test not available: {}".format(msg))
            return
        manager = cpd.PhononDOS_DB(db_name)
        omega = np.linspace(0.0,100.0,10)
        dos = omega**2
        atID = 1
        name = "firstEntry"
        id = 0
        manager.save( name=name, atID=atID, omega_e=omega, dos_e=dos )

        # Load back in based on id
        res = manager.get(id=id)
        self.assertTrue( res["id"] == id )
        self.assertTrue( res["atID"] == atID )
        self.assertTrue( np.array_equal(res["omega_e"], omega) )
        self.assertTrue( np.array_equal(res["dos_e"],dos) )

        # Load back in based on atomID
        res = manager.get(atID=atID)
        self.assertTrue( res["id"] == id )
        self.assertTrue( res["atID"] == atID )
        self.assertTrue( np.array_equal(res["omega_e"], omega) )
        self.assertTrue( np.array_equal(res["dos_e"],dos) )

        # Load back in based on name
        res = manager.get(name=name)
        self.assertTrue( res["id"] == id )
        self.assertTrue( res["atID"] == atID )
        self.assertTrue( np.array_equal(res["omega_e"], omega) )
        self.assertTrue( np.array_equal(res["dos_e"],dos) )

    def __del__(self):
        # Clean up
        os.remove(db_name)
