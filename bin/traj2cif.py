#!/usr/bin/env python
import sys
import os
from ase.io import Trajectory, write

def main(fname=None, ext="cif"):
    if fname is None:
        raise ValueError("No file name given")

    folder_name = fname.rpartition(".")
    os.mkdir(folder_name)
    traj = Trajectory(fname)
    for i,atom in enumerate(traj):
        atoms_file = folder_name+"/atoms{}.{}".format(i,ext)
    print ("Atoms written to {}".format(folder_name))

def print_usage(known_args):
    print ("Usage:")
    kw_string = ""
    for arg in known_args:
        kw_string += "--{}=<...> ".format(arg)
    print ("traj2cif.py {}".format(kw_string))

if __name__ == "__main__":
    known_args = ["fname","ext"]
    kwargs = {}
    for arg in sys.argv[1:]:
        arg_found = False
        for cand_arg in known_args:
            if arg.find("--{}=".format(cand_arg)) != -1:
                kwargs[cand_arg] = arg.split("--{}=".format(cand_arg))[1]
                arg_found = True
                break
        if not arg_found:
            def print_usage(known_args)
            raise ValueError("Unknown argument {}. Known arguments: {}".format(arg,known_args))

    main(**kwargs)
