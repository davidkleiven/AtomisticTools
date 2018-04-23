import numpy as np
from matplotlib import pyplot as plt
from ase.db import connect
from ase.ce import Evaluate

class CVScoreHistory(object):
    def __init__( self, setting=None, penalization="L1", select_cond=None ):
        self.setting = setting
        self.penalization = penalization
        self.select_cond = select_cond
        self.result = None

    def get_generations( self ):
        """
        Compute the CV score and RMSE error
        """
        db = connect( self.setting.db_name )
        select_cond = [("converged","=","1")]
        if ( self.select_cond is not None ):
            select_cond += self.select_cond

        all_generations = []
        for row in db.select(select_cond):
            try:
                gen = row["gen"]
            except KeyError as exc:
                print (str(exc))
                continue
            if ( gen not in all_generations ):
                all_generations.append( gen )

        all_generations.sort()
        print ("Generations: {}".format(all_generations))
        del all_generations[0]
        return all_generations

    def reset(self):
        """
        Reset such that the history is recomputed
        """
        self.result = None

    def get_history( self, lambdas=None ):
        """
        Compute the history
        """

        if ( lambdas is None ):
            raise ValueError( "No lambdas given!" )

        try:
            N = len(lambdas)
        except:
            raise TypeError( "Lambdas has to be a list" )

        ecis = {}
        gens = self.get_generations()

        cv_gen = []
        rmse_gen = []
        num_structs = []
        for i,gen in enumerate(gens):
            scond = [("gen","<=",gen)]
            print ("Current generation: {} ({}%)".format(gen,int(100*i/len(gens)) ) )
            if ( self.select_cond is not None ):
                scond += self.select_cond
            cvs = []
            for lamb in lambdas:
                evaluator = Evaluate( self.setting, lamb=float(lamb), penalty=self.penalization, select_cond=scond )
                cvs.append( evaluator._cv_loo() )
            indx = np.argmin(cvs)
            print ("Selected penalization value. Indx: {}. Value: {}".format(indx,lambdas[indx]))
            evaluator = Evaluate( self.setting, lamb=float(lambdas[indx]), penalty=self.penalization, select_cond=scond )
            current_ecis = evaluator.get_cluster_name_eci_dict
            small_cv = evaluator._cv_loo()
            evaluator._get_dft_energy_per_atom()
            evaluator._get_e_predict()
            rmse = evaluator.rmse()
            #rmse = small_cv
            cv_gen.append( small_cv )
            rmse_gen.append( rmse )
            num_structs.append( evaluator.cf_matrix.shape[0] )
            for key,value in current_ecis.iteritems():
                if ( key not in ecis.keys() ):
                    ecis[key] = {"gen":[],"value":[]}
                ecis[key]["gen"].append(gen)
                ecis[key]["value"].append(value)

        self.result = {
            "cv":cv_gen,
            "rmse":rmse_gen,
            "gen":gens,
            "ecis":ecis,
            "num_structs":num_structs
        }
        return self.result

    def plot(self, max_eci_per_plot=10 ):
        """
        Plots the result
        """
        if ( self.result is None ):
            raise ValueError( "The history has not been computed yet. Call get_history() first" )

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot( self.result["gen"], np.array(self.result["cv"])*1000.0, marker="o", mfc="none", label="CV" )
        ax.plot( self.result["gen"], np.array(self.result["rmse"])*1000.0, marker="x", label="RMSE" )
        ax.set_xlabel( "Generation" )
        ax.set_ylabel( "CV/RMSE (meV/atom)" )
        ax.legend(loc="best",frameon=False)
        ax2 = ax.twinx()
        ax2.plot( self.result["gen"], self.result["num_structs"], "--", ls="steps", color="#bdbdbd")
        ax2.set_ylabel( "Number of structures" )

        # Plot how the ECIs evolve
        n_ecis = len(self.result["ecis"].keys() )
        n_plots = int( n_ecis/max_eci_per_plot )
        n_cols = int( np.sqrt(n_plots) ) + 1

        if ( n_plots == 1 ):
            n_rows = 1
        else:
            n_rows = int(n_plots/(n_cols-1) )

        fig_eci, ax_eci = plt.subplots( nrows=n_rows, ncols=n_cols, sharex=True )

        counter = 0
        plot_no = 0
        for key,value in self.result["ecis"].iteritems():
            plot_no = int(counter/max_eci_per_plot)
            col = plot_no%n_cols
            row = int(plot_no/n_cols)
            ax_eci[row,col].plot( value["gen"], value["value"], label=key )
            counter += 1
            if ( counter%max_eci_per_plot == 0 ):
                plot_no += 1
        return fig, fig_eci
