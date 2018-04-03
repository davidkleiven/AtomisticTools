import unittest
from atomtools.eos import RedlichKister

class TestRedlichKister(object):
    def test_nothrow(object):
        no_throw = True
        msg = ""
        try:
            mg_conc = np.linspace(0.0,1.0,20)
            al_conc = 1.0-mg_conc
            comp = {
                "Al":al_conc,
                "Mg":mg_conc
            }
            free_energy = al_conc*mg_conc
            ref_eng = {
                "Al":-1.0,
                "Mg":-0.5
            }
            redkist = RedlichKister( ref_eng )
            redkist.fit( free_energy, comp )
            excess = redkist.eval( comp )
        except Exception as exc:
            msg = str(exc)
            no_throw = False
        self.assertTrue( no_throw, msg=msg )

if __name__ == "__main__":
    unittest.main()
