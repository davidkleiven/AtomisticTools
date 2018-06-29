from setuptools import setup, find_packages

setup(
    name="atomtools",
    author="David Kleiven",
    version=1.0,
    packages=find_packages(),
    scripts=["bin/viewdb.py","bin/traj2cif.py", "bin/center_atom.py"]
)
