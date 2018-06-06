from ase.db import connect
from ase.optimize import BFGS
from ase.constraints import StrainFilter
from ase.optimize.precon import PreconLBFGS
from atomtools.ase import align_direction_with_z
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize

class ConstituentStrain(object):
    """
    Class for computing the constituent strain according to
    Mixed-basis cluster expansion for thermodynamics of bcc alloys, Volker Blum and Alex Zunger, PHYSICAL REVIEW B 70, 2004
    """
    def __init__( self, atoms=None, db_name=None ):
        self.atoms = atoms
        self.orig_atoms = self.atoms.copy()
        self.orig_atoms.set_calculator( self.atoms._calc )
        self.db_name = db_name

    def run_one_in_plane_distance( self, cell_length_x=2.0, cell_length_y=2.0, direction=(0,0,1), smax=0.003 ):
        # Rotate the atoms object such that the z-direction points along target_z_direction
        self.atoms = align_direction_with_z( self.atoms, direction=direction )
        lengths_angles = self.atoms.get_cell_lengths_and_angles()
        ratio_x = cell_length_x/lengths_angles[0]
        ratio_y = cell_length_y/lengths_angles[1]
        cell = self.atoms.get_cell()
        cell[0,:] *= ratio_x
        cell[1,:] *= ratio_y
        self.atoms.set_cell( cell.T, scale_atoms=True )

        strfilter = StrainFilter( self.atoms, mask=[0,0,1,1,1,1] )
        relaxer = BFGS(self.atoms)
        V = self.atoms.get_volume()
        relaxer.run( fmax=smax*V )

        # Store the results to ase db
        db = connect( self.db_name )
        kvp = {
            "cell_length_x":cell_length_x,
            "cell_length_y":cell_length_y,
            "direction_x":direction[0],
            "direction_y":direction[1],
            "direction_z":direction[2]
        }
        db.write( self.atoms,key_value_pairs=kvp )

    def full_relaxation( self, fmax=0.025, smax=0.003 ):
        """
        Perform a full relaxation of the system
        """
        if ( len(self.atoms) == 1 ):
            strfilter = StrainFilter(self.atoms)
            relaxer = BFGS(strfilter)
            fmax = smax*self.atoms.get_volume()
            relaxer.run(fmax=fmax)
        else:
            relaxer = PreconLBFGS(self.atoms, variable_cell=True )
            relaxer.run( fmax=fmax, smax=smax )
        db = connect( self.db_name )
        db.write( self.atoms, key_value_pairs={"full_relaxation":True} )

    def run( self, cell_length_x=None, cell_length_y=None, direction=(0,0,1), smax=0.003 ):
        """
        Run all cell lengths
        """
        if ( cell_length_y is None ):
            cell_length_y = cell_length_x
        if ( cell_length_x is None ):
            raise ValueError( "No cell length is given!" )

        for ax,ay in zip(cell_length_x,cell_length_y):
            try:
                self.run_one_in_plane_distance( cell_length_x=ax, cell_length_y=ay, direction=direction, smax=smax )
            except Exception as exc:
                print (str(exc))
                print ("Proceeding to the next")
            self.atoms = self.orig_atoms.copy()
            self.atoms.set_calculator( self.orig_atoms._calc )

    def get_energy_in_plane_distance( self, direction=None ):
        db = connect( self.db_name )
        energies = {}
        a_in_plane = {}

        for row in db.select( direction_x=direction[0], direction_y=direction[1], direction_z=direction[2] ):
            if ( row.formula not in energies.keys() ):
                energies[row.formula] = []
                a_in_plane[row.formula] = []
            energies[row.formula].append( row.energy/row.natoms )
            a_in_plane[row.formula].append( row.cell_length_x )
        return energies,a_in_plane
    def plot_strain_energies( self, direction=None, ax=None ):
        if ( ax is None ):
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)

        energies, a_in_plane = self.get_energy_in_plane_distance( direction=direction )

        for key in energies.keys():
            E = energies[key]
            a = a_in_plane[key]
            srt_indx = np.argsort(a)
            print (srt_indx)
            a = np.array( [a[indx] for indx in srt_indx] )
            E = np.array( [E[indx] for indx in srt_indx] )
            ax.plot( a, E-E[0], marker="o", label=key )
        ax.legend( frameon=False, loc="best" )
        return fig

    def coherency_strain_energy( self, direction=None ):
        """
        Computes the coherency strain energy
        """
        energies, a_in_plane = self.get_energy_in_plane_distance( direction=direction )
        if ( len(energies.keys()) > 2 ):
            msg = "Systems with more than two elements are not supported!\n"
            msg += "Elements found in database {}".format(energies.keys())
            raise ValueError( msg )

        key1 = energies.keys()[0]
        key2 = energies.keys()[1]
        i1 = interp1d( a_in_plane[key1], energies[key1], bounds_error=False, fill_value="extrapolate", kind="linear" )
        i2 = interp1d( a_in_plane[key2], energies[key2], bounds_error=False, fill_value="extrapolate", kind="linear" )

        # Find the minimum of each function
        def cost_i1(a):
            return i1(a)
        def cost_i2(a):
            return i2(a)

        x01 = np.argmin(energies[key1])
        x02 = np.argmin(energies[key2])
        res1 = minimize(cost_i1,x0=x01)
        res2 = minimize(cost_i2,x0=x02)
        ref1 = i1(res1["x"])
        ref2 = i2(res2["x"])
        i1 = interp1d( a_in_plane[key1], np.array(energies[key1])-ref1, bounds_error=False, fill_value="extrapolate", kind="linear" )
        i2 = interp1d( a_in_plane[key2], np.array(energies[key2])-ref2, bounds_error=False, fill_value="extrapolate", kind="linear" )

        def cost_func( a, conc ):
            return i1(a)*conc + (1.0-conc)*i2(a)

        concs = np.linspace(0.0,1.0,20)
        a_in_plane_min = []
        E_CS = []
        for conc in concs:
            x0 = conc*x01 + (1.0-conc)*x02
            res = minimize( cost_func, x0=x0, args=(conc,) )
            a_in_plane_min.append(res["x"])
            E_CS.append( res["fun"] )
        return concs, E_CS
