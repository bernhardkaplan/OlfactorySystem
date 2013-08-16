import re
import os
import sys
import numpy as np
import pylab
import simulation_parameters
from FigureCreator import plot_params
import json


class AnalyseResponseCurves(object):
    """
    Calculates for a set of response curves the use defined fitness
    """

    def __init__(self, params):
        self.params = params


    def sweep_orn_response_fitness(self):
        match_files = 'orn_response_curve_(\d+).dat'
#        match_files = 'orn_response_curve_%d.dat'
        file_dict = self.find_files(match_files)
        n_matches = len(file_dict.keys())
        assert (n_matches > 0), 'No files found matching the pattern %s in folder %s\nCorrect folder in simulation_parameters.py?' % (match_files, self.params['other_folder'])
        fitness_values = np.zeros(n_matches)
        print 'n_matches', n_matches
        for i in xrange(n_matches):
            fn = file_dict[i]
            d = np.loadtxt(fn)
            fnv = self.get_orn_response_fitness(d)
            print 'Analysing', fn, fnv
            fitness_values[i] = fnv

        n_best = 30
        best_sims = fitness_values.argsort()[-n_best:]
        print 'Best simulations and fitness values'
        print best_sims
        print fitness_values[best_sims]

        fig_fns = []
        display_command = 'ristretto '
        for i_ in xrange(n_best):
            fig_fn = self.params['figure_folder'] + '/hand_tuned_%d.png' % best_sims[i_]
            fig_fns.append(fig_fn)
            display_command += '%s ' % fig_fn
            print fig_fn, 
        print ' '
        output_fn = self.params['other_folder'] + '/good_figure_filenames.json' 
        output_f = file(output_fn, 'w')

#        print 'debug', type(fig_fns)
#        print 'Writing good figure filenames to:', output_fn
        os.system(display_command)
#        json.dump(output_f, fig_fns)


    def get_orn_response_fitness(self, xy_data):
        """
        Returns the fitness value for a set of response curves.

        Keyword arguments:
        xy_data -- numpy.ndarray
        xy_data[:, 0] = concentration values
        xy_data[:, i] = output-rate of curve i
        """
        fitness = 0.
        
        x_data = xy_data[:, 0]
        y_data = xy_data[:, 1:]
        n_x = x_data.size
        n_curves = xy_data[0, 1:].size # the number of curves to be analysed

        # define some fixed characteristica to be applied to the set of curves
        # for each curve there are some individual criteria, i.e. how the curve should look like
        y_min, y_max = .05, .9 # percentage of the curve's global y-min (max) which marks the x_min, x_max points
        x_mins, x_maxs = (-1) * np.ones(n_curves), (-1) * np.ones(n_curves) # negative to check if they've already been set
        y_max_global = np.zeros(n_curves)

        # get the curve characteristics
        for curve_ in xrange(n_curves):
            y_max_curve = y_max * y_data[:, curve_].max()
            y_min_curve = y_min * y_max_curve
            y_max_global[curve_] = y_data[:, curve_].max()
            for i_ in xrange(n_x):
                x, y = x_data[i_], y_data[i_, curve_]
                if y > y_min_curve and (x_mins[curve_] == -1.):
                    x_mins[curve_] = x
                if y > y_max_curve and (x_maxs[curve_] == -1.):
                    x_maxs[curve_] = x


        punishment = 20
        # compare the curves within the set
        for curve_ in xrange(1, n_curves):
            # compare global maxima
            dy = y_max_global[curve_] - y_max_global[curve_-1]
            if dy < 0: # quite bad
#                fitness += dy
                fitness -= punishment
            else:
                fitness += punishment
            # compare the point when curves cross the y_min threshold
            dx_mins = x_mins[curve_] - x_mins[curve_-1]
            if dx_mins < 0:
                fitness -= punishment
            else:
                fitness += punishment
            # compare the point when curves cross the y_max threshold
            dx_maxs = x_maxs[curve_] - x_maxs[curve_-1]
            if dx_maxs < 0:
                fitness -= punishment
            else:
                fitness += punishment
#            print 'curve_ xmin, xmax, dy', curve_, x_mins[curve_], x_maxs[curve_], dy

        # compare the highest and smallest ymax values of all curves
        dy_max_global = np.max(y_max_global) - np.min(y_max_global)
        fitness -= dy_max_global / (np.min(y_max_global) + 1e-6)
        # punish large differences in maximum output rates
        if dy_max_global > 20:
            fitness -= punishment
        # punish too small output rates
        if np.mean(y_max_global) < 20:
            fitness -= punishment
        # punish too high max output rates
        if np.max(y_max_global) > 160:
            fitness -= punishment


        # compare the points when curves cross their y_min
        dx_min_global = x_mins[0] - x_mins[-1]
        fitness += dx_min_global * 10.

        return fitness




    def find_files(self, pattern, folder=None):
        """
        Returns a dictionary with all absolute file names 
        which match a certain pattern in a given folder.
        The key of the file dictionary is the simulation ID 
        and the value is the absolute path.
        """

        if folder == None:
            folder = self.params['other_folder']
        file_dict = {}
        for f in os.listdir(folder):
            m = re.match(pattern, f)
            if m:
                sim_id = int(m.groups()[0])
                file_dict[sim_id] = os.path.abspath(folder + '/' + f)
        return file_dict


    def create_new_parameters_from_selected_sets(self, list_of_sim_ids, param_log_fn_base=None):
        """
        list_of_sim_ids: list of integer values with the simulation ids that contain good parameter sets
        """
        if param_log_fn_base == None:
            param_log_fn_base = self.params['params_folder'] + '/params_handtuning_'

        for i_, sim_id in enumerate(list_of_sim_ids):
            fn = param_log_fn_base + '%d.txt' % sim_id
            print 'debug', os.path.exists(fn)




if __name__ == '__main__':

    param_tool = simulation_parameters.parameter_storage(use_abspath=True)
    #param_tool = simulation_parameters.parameter_storage(use_abspath=False)
    params = param_tool.params

    ARC = AnalyseResponseCurves(params) 
    ARC.sweep_orn_response_fitness()
