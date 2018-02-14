import numpy as np
import matplotlib as mpl
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["axes.unicode_minus"] = False
#mpl.rcParams["font.size"] = 18
from matplotlib import pyplot as plt
import copy

class CovariancePlot( object ):
    def __init__( self, evaluator, constant_term_column=None ):
        self.eval = evaluator
        self.cf_mat = evaluator.cf_matrix
        self.constant_term_column = constant_term_column

    def plot( self ):
        indx = []
        cf_mat = copy.deepcopy( self.cf_mat )
        if ( not self.constant_term_column is None ):
            cf_mat = np.delete( cf_mat, self.constant_term_column, axis=1 )
        for i in range(0,self.cf_mat.shape[0]):
            if ( np.std(self.cf_mat[i,:]) == 0 ):
                indx.append(i)
        cf_mat = np.delete( cf_mat, indx, axis=0 )
        correlation_matrix = np.corrcoef( cf_mat )

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.imshow( correlation_matrix, interpolation="none", cmap="coolwarm")
        cbar = fig.colorbar(im)
        cbar.set_label( "Structure correlation" )
        ax.set_xlabel( "Structure number" )
        ax.set_ylabel( "Structure number" )
        return fig

    def plot_corr_func_coverage(self,nbins=50):
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        all_data = {}
        labels = {
            1:"Singlets",
            2:"Doublets",
            3:"Triplets",
            4:"Quads",
        }
        prec = np.linalg.inv(self.cf_mat.T.dot(self.cf_mat))
        prec_num = np.trace(prec)/self.cf_mat.shape[0]
        print (prec_num)
        for i in range(self.cf_mat.shape[1]):
            cname = self.eval.cluster_names[i]
            size = int(cname[1])
            if ( not size in all_data.keys() ):
                all_data[size] = []
            all_data[size] += list(self.cf_mat[:,i])
        bottom = np.zeros(nbins)
        bins = np.linspace(-1.0,1.0,nbins+1)
        bins = bins
        delta = bins[1]-bins[0]
        for key,value in all_data.iteritems():
            if ( key==0 ):
                continue
            hist,bins = np.histogram( value, bins=bins )
            hist = hist.astype(np.float64)
            hist /= float(np.sum(np.abs(hist) ))
            ax.bar( bins[1:], hist, 0.8*delta, bottom=bottom, label=labels[key] )
            bottom += hist
        ax.legend( loc="best", frameon=False )
        ax.set_xlabel( "Correlation function" )
        ax.set_ylabel( "Relative occurence" )
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        return fig
