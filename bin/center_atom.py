#!/usr/bin/env python

"""
This script translate an atom to the center of the cell,
then wraps all the atoms inside the unitcell
"""
import sys
from ase.io import read, write

def main(argv):
    allowed_args = ["--fname=<filename> --atom=<indx of atoms to put at center>"]

    for arg in argv:
        if arg.find("--fname=") != -1:
            fname = arg.split("--fname=")[1]
        elif arg.find("--atom=") != -1:
            atom_indx = int(arg.split("--atom=")[1])
        else:
            print("Unknown argument!")
            print("Has to be one of {}".format(allowed_args))
            return

    ext = fname.rpartition(".")[-1]
    out_fname = fname.rpartition(".")[0]+"_{}.{}".format(atom_indx, ext)

    atoms = read(fname)
    cell = atoms.get_cell()
    diag = 0.5*(cell[0,:] + cell[1,:] + cell[2,:])
    atoms.translate(diag-atoms[atom_indx].position)
    atoms.wrap()
    write(out_fname, atoms)
    print("Centered atoms object written to {}".format(out_fname))

if __name__ == "__main__":
    main(sys.argv[1:])
