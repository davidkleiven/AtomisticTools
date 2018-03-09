from matplotlib import pyplot as plt
import numpy as np
from ase.units import kB, kg, _hbar, _c, eV, _amu, J, Angstrom
from scipy.optimize import minimize
from ase.data import atomic_masses, atomic_numbers

class EquationOfState(object):
    def __init__( self, volume, energy, debye_scheme="mjs" ):
        self.volume = volume
        self.energy = energy
        self.tot_mass = None
        self.natoms = 1
        allowed_debye_schemes = ["none","mjs"]
        if ( debye_scheme not in allowed_debye_schemes ):
            raise ValueError( "Debye Scheme has to be one of {}".format(allowed_debye_schemes) )
        self.debye_scheme = debye_scheme

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
        if ( self.tot_mass is None ):
            raise ValueError( "Average mass is not known" )
        rho = self.tot_mass/V
        return rho

    def density_g_per_cm3( self, V ):
        """
        Returns the density in g/cm^3
        """
        rho = self.density(V)
        rho /= kg
        rho *= 1E30
        return rho/1000.0

    def set_average_mass( self, atoms ):
        """
        Computes the density given a number dictionary of atoms
        """
        n_tot = 0
        for key,value in atoms.iteritems():
            n_tot += value
        self.natoms = n_tot
        tot_mass = 0.0
        for key,value in atoms.iteritems():
            tot_mass += value*atomic_masses[atomic_numbers[key]]
        self.tot_mass = tot_mass
        #print (self.avg_mass)

    def debye_frequency( self, V ):
        """
        Compute the Debye Frequency in rad/s of a given assuming isotropic speed of sound
        """
        B = self.bulk_modulus(V)
        rho = self.density(V)
        rho /= kg # Unit of rho is now kg/(angstrom^3)
        B /= J # J/angstrom^3
        omega_D = (6.0*np.pi**2)**(1.0/3.0) * V**(-1.0/3.0) * np.sqrt( B/rho )*1E10 # Unit of omega_D: rad/s
        v_sound = self.speed_of_sound(V)
        number_density = self.natoms/V
        omega_D = (6.0*np.pi**2 *number_density)**(1.0/3.0) * v_sound*1E10
        omega_D *= _hbar*J

        if ( self.debye_scheme == "none" ):
            return omega_D # In eV
        elif ( self.debye_scheme == "mjs" ):
            # See: Moruzzi, V.; Janak, J. & Schwarz, K. Calculated thermal properties of metals Physical Review B, APS, 1988, 37, 790
            return 0.617*omega_D
        raise ValueError( "Unknown debye scheme!" )

    def debye_temperature( self, V ):
        """
        Computes the Debye temperature
        """
        return self.debye_frequency(V)/kB

    def speed_of_sound( self, V ):
        """
        Computes the speed of sound in meter per second
        """
        B = self.bulk_modulus(V)
        rho = self.density(V)
        rho /= kg
        B /= J
        return np.sqrt(B/rho)

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
        return res["fun"], res["x"][0]

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

        Emin = self.minimum_energy()[0]/natoms
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
