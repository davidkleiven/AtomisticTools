#!/usr/bin/env python
import glob
import sys
from ase.io import write
import gpaw as gp


def main(folder):
    for fname in glob.glob(folder+"/*.gpw"):
        atoms, calc = gp.restart(fname)
        out_fname = fname.rsplit(".")[0]+".xyz"
        write(out_fname, atoms)
        print("Saved atoms to {}".format(out_fname))


if __name__ == "__main__":
    folder = sys.argv[1]
    main(folder)
