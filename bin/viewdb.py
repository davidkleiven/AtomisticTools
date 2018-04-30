#!/usr/bin/env python
import sys
from ase.db import connect
from ase.visualize import view

def show(db_name,uid):
    db = connect( db_name )
    atoms = db.get_atoms(id=uid)
    view(atoms)

def main(argv):
    for arg in argv:
        if ( arg.find("help") != -1 or arg.find("-h") != -1 ):
            print ("Usage: viewdb.py <db_name> <uid>")
            return
    db_name = argv[0]
    uid = int(argv[1])
    show(db_name,uid)

if __name__ == "__main__":
    main(sys.argv[1:])
