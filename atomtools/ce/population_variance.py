from cemc import CE
from ase.ce.corrFunc import CorrFunction
import numpy as np
import time
import json
from matplotlib import pyplot as plt
import h5py as h5

class PopulationVariance( object ):
    """
    Object that estimates the covariance and the mean of all the
    correlation functions in the population
    """
    def __init__( self, BC ):
        self.bc = BC
        cf_obj = CorrFunction(self.bc)
        self.init_cf = cf_obj.get_cf( self.bc.atoms )
        ecis = {key:1.0 for key in self.init_cf.keys()} # ECIs do not matter here
        self.calc = CE( self.bc, ecis, initial_cf=self.init_cf )
        self.bc.atoms.set_calculator(self.calc)
        self.elements = self.bc.basis_elements[0] # This should be updated to handle different site types
        self.status_every_sec = 30
        N = len(ecis.keys())-1
        self.cov_matrix = np.zeros((N,N))
        self.exp_value = np.zeros(N)

    def swap_random_atoms( self ):
        """
        Changes the symbol of a random atom
        """
        indx = np.random.randint(low=0,high=len(self.bc.atoms))
        symb = self.bc.atoms[indx].symbol
        new_symb = self.elements[np.random.randint(low=0,high=len(self.elements))]
        system_change = [(indx,symb,new_symb)]
        self.bc.atoms._calc.calculate( self.bc.atoms, ["energy"], system_change )

    def estimate( self, n_probe_structures=10000, fname="" ):
        """
        Estimate the covariance matrix
        """

        if ( fname != "" ):
            if ( not fname.endswith(".json") ):
                raise ValueError( "The program is going to write a json file so the file extension should be .json" )

        step = 0
        cfs = []
        cur_time = time.time()
        while ( step < n_probe_structures ):
            step += 1
            if ( time.time()-cur_time > self.status_every_sec ):
                print ("Step {} of {}".format(step,n_probe_structures) )
                cur_time = time.time()
            for i in range(len(self.bc.atoms)):
                self.swap_random_atoms()
            new_cfs = self.calc.get_cf()
            cf_array = [new_cfs[key] for key in self.init_cf.keys() if key != "c0"] # Make sure that the order is the same as in init_cf
            self.update_cov_matrix( cf_array )

        cfs = np.array( cfs )
        cov = self.cov_matrix/n_probe_structures
        mu = self.exp_value/n_probe_structures
        cov -= np.outer(mu,mu)
        return cov,mu

    def array2dict( self, cov, mu ):
        """
        Converts an array to a dictionary
        """
        # Create dictionaries
        keys = self.init_cf.keys()
        del keys[keys.index("c0")]
        mu_dict = {keys[i]:mu[i] for i in range(len(keys))}

        cov_dict = {key:{} for key in keys}
        for i in range(len(keys)):
            for j in range(len(keys)):
                cov_dict[keys[i]][keys[j]] = cov[i,j]
        return cov_dict, mu_dict


    def update_cov_matrix( self, new_cfs ):
        """
        Updates the covariance matrix
        """
        outer = np.outer( new_cfs, new_cfs )
        self.cov_matrix += outer
        self.exp_value += np.array( new_cfs )

    def save( self, eigval, eigvec, fname="eigenvectors.h5", fraction=0.95 ):
        """
        Store the lowest fraction of the eigenvectors
        """
        srt_indx = np.argsort( eigval )[::-1]
        eigval_srt = [eigval[indx] for indx in srt_indx]
        eigvec_srt = np.array( [eigvec[:,indx] for indx in srt_indx] ).T
        cumsum_eig = np.cumsum( eigval_srt )
        tot_sum = np.sum(eigval_srt)

        relative_sum = cumsum_eig/tot_sum

        max_indx = np.argmin( np.abs(relative_sum-fraction) )
        matrix = eigvec_srt[:,:max_indx]
        eigen = eigval_srt[:max_indx]
        with h5.File( fname, 'w' ) as hf:
            dset_vec = hf.create_dataset( "eigenvectors", data=matrix )
            dset_eigen = hf.create_dataset( "eigenvalues", data=eigen )
            dset_eigen.attrs["fraction"] = fraction

        print ( "Result written to {}".format(fname) )

    def diagonalize( self, cov, plot=False ):
        """
        Diagonalize the covariance matrix
        """
        eigval, eigvec = np.linalg.eigh( cov )

        # Sort accorting to eigenvalues
        srt_indx = np.argsort( eigval )[::-1]
        eigval_srt = [eigval[indx] for indx in srt_indx]
        eigvec_srt = np.array( [eigvec[:,indx] for indx in srt_indx] ).T

        cumsum_eig = np.cumsum( eigval_srt )
        tot_sum = np.sum(eigval_srt)

        if ( plot ):
            x_val = np.arange(len(eigval))
            fig = plt.figure()
            ax = fig.add_subplot(1,1,1)
            ax.plot( x_val, cumsum_eig/tot_sum, ls="steps" )
            ax.set_xlabel( "Number of eigenvectors" )
            ax.set_ylabel( "Normalized variance" )
        return eigval, eigvec

    def plot_eigen( self, eigenvalues, eigenvectors ):
        """
        Creates a visualization of the eigenvectors and the eigenvalues
        NOTE: Eigenvalues and eigenvectors must NOT be sorted. They have
              to be given in the same order as returned by diagonalize
        """
        keys = self.init_cf.keys()
        del keys[keys.index("c0")]
        srt_indx = np.argsort(eigenvalues)[::-1]
        key_srt = [keys[indx] for indx in srt_indx]
        eigval_srt = [eigenvalues[indx] for indx in srt_indx]
        eigvec_srt = np.array( [eigenvectors[:,indx] for indx in srt_indx] ).T

        grid_kw = {"hspace":0.0,"height_ratios":[1,4]}
        fig, ax = plt.subplots(nrows=2,sharex=True,gridspec_kw=grid_kw)
        x = np.arange(len(eigvec_srt))
        cumsum = np.cumsum(eigval_srt)
        tot_sum = np.sum(eigval_srt)
        ax[0].plot( x, cumsum/tot_sum, ls="steps")
        keys_srt, eigvec_srt = self.sort_eigenvectors_by_size( keys, eigvec_srt )
        #ax[1].imshow( eigvec_srt, cmap="coolwarm", aspect="auto" )
        print (eigvec_srt[:,0])

        # Separation lines to separate different sizes
        hlines = []
        size = 2
        for i,name in enumerate(keys_srt):
            prefix = "c{}".format(size)
            if ( name.startswith(prefix) ):
                hlines.append(i)
                size += 1

        for line in hlines:
            ax[1].axhline( line )
        return fig

    def sort_eigenvectors_by_size( self, keys, eigvec ):
        """
        Sort the eigenvector values
        """
        srt_indx = np.argsort(keys)
        keys_srt = [keys[indx] for indx in srt_indx]
        eigvec_srt = np.array( [eigvec[indx,:] for indx in srt_indx] )
        return keys_srt, eigvec_srt
