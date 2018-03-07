from matplotlib import pyplot as plt
import numpy as np
from ase.units import kB, _hbar, _c, eV
from scipy.optimize import minimize
from ase.data import atomic_masses, atomic_numbers

class EquationOfState(object):
    def __init__( self, volume, energy ):
        self.volume = volume
        self.energy = energy
        self.avg_mass = None

    def evaluate( self, V ):
        raise NotImplementedError( "This function has to be implemented in subclasses" )

    def deriv( self, V ):
        raise NotImplementedError( "Derivatives has to be implemented in subclasses" )

    def double_deriv( self, V ):
        raise NotImplementedError( "Double derivative has to be implemented in subclasses" )

    def plot( self, latex=False ):
        """
        Plots the result
        """
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        V = np.linspace( 0.95*np.min(self.volume), 1.05*np.max(self.volume), 100 )
        fitted = self.evaluate(V)
        ax.plot( V, fitted )
        ax.plot( self.volume, self.energy, "o", mfc="none" )
        if ( latex ):
            ax.set_xlabel( "Volume (\$\SI{}{\\angstrom^3}\$)")
        else:
            ax.set_xlabel( r"Volume (\AA^3)" )
        ax.set_ylabel( "Energy (eV)" )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        return fig

    def bulk_modulus( self, V ):
        """
        Computes the bulk modulus as a function of volume
        """
        B = self.double_deriv(V)*V
        if ( isinstance(B,np.ndarray) ):
            B[B<0.0] = 0.0
        return B

    def density( self, V ):
        """
        Compute the density based on the reference density
        """
        if ( self.avg_mass is None ):
            raise ValueError( "Average mass is not known" )
        rho = self.avg_mass/V
        return rho

    def set_average_mass( self, atoms ):
        """
        Computes the density given a number dictionary of atoms
        """
        n_tot = 0
        for key,value in atoms.iteritems():
            n_tot += value
        avg_mass = 0.0
        for key,value in atoms.iteritems():
            avg_mass += value*atomic_masses[atomic_numbers[key]]
        self.avg_mass = avg_mass/n_tot

    def debye_frequency( self, V ):
        """
        Compute the Debye Frequency in rad/s of a given assuming isotropic speed of sound
        """
        B = self.bulk_modulus(V)
        rho = self.density(V)
        u =  931.4941E6 # eV/c^2
        rho *= u
        omega_D = (6.0*np.pi**2)**(1.0/3.0) * V**(-1.0/3.0) * np.sqrt( B/rho )
        omega_D *= _c*1E10*_hbar
        eV_si = 1.6021766208E-19 # eV/J
        omega_D /= eV_si
        return omega_D # In eV

    def phonon_free_energy_high_temp( self, debye_freq, T ):
        """
        Computes the phonon free energy in the high temperature limit
        """
        return kB*T*( 3.0*np.log(debye_freq/(kB*T)) - 1.0 )

    def minimum_energy( self ):
        """
        Compute the minimum energy
        """
        indx = np.argmin(self.energy)
        V0 = self.volume[indx]
        res = minimize( self.evaluate, V0 )
        return res["fun"]

    def beta_elastic_vib_free_energy( self, T, vol_curve=None, natoms=1 ):
        """
        Compute the total free energy elastic + vibration
        """
        if ( vol_curve is None ):
            vol_curve = self.volume_temperature( T, natoms )

        fvib = []
        elastic = self.evaluate( vol_curve )/natoms
        for V,temp in zip(vol_curve,T):
            debye = self.debye_frequency(V)
            fvib.append( self.phonon_free_energy_high_temp(debye,temp) )

        Emin = self.minimum_energy()/natoms
        E = elastic-Emin + np.array(fvib)
        if ( len(T) == 1 ):
            return E[0]/(kB*T[0])
        return E/(kB*np.array(T))

    def volume_temperature( self, T, natoms ):
        """
        Computes the volume as a function of temperature by minimizing the
        Free Energy of elastic + vibration
        """
        min_indx = np.argmin(self.energy)
        V0 = self.volume[min_indx]
        volumes = []
        for temp in T:
            res = minimize( minization_vol_temp_curve, V0, args=(self,natoms,temp) )
            volumes.append(res["x"][0])
            V0 = res["x"]
        return np.array(volumes)

    def linear_thermal_expansion_coefficient( self, T, natoms=1, vol_curve=None ):
        """
        Computes the linear thermal expansion coefficient.

        Relation between linear (alpha_L) and volume (alpha_V) expansion
        alpha_V = 3*alpha_L
        """
        if ( vol_curve is None ):
            # Volume curve not given so we have to compute it
            vol_curve = self.volume_temperature( T, natoms )

        X = np.zeros((len(T),3))
        X[:,0] = 1.0
        X[:,1] = T
        X[:,2] = T**2
        res, residuals,rank,s = np.linalg.lstsq( X, vol_curve )
        alpha_V = (res[1]+2.0*res[2]*T)/vol_curve
        alpha_L = alpha_V/3.0
        return alpha_L

def minization_vol_temp_curve( V, eos, natoms, temperature ):
    elastic_energy = eos.evaluate(V)/natoms
    debye_freq = eos.debye_frequency( V )
    F_vib = eos.phonon_free_energy_high_temp( debye_freq, temperature )
    return elastic_energy + F_vib
