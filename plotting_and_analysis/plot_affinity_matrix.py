import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import pylab
import numpy as np
import sys
import os
import simulation_parameters

class OdorantAffininityPlotter(object):
    def __init__(self, params):
        self.params = params

    def plot_affinity_matrix(self, path=None):
        if path == None:
            path = self.params['activation_matrix_fn']
            print 'Loading', path
        d = np.loadtxt(path)
        self.activation_matrix = d


        plot_params = {'backend': 'png',
                      'axes.labelsize': 20,
                      'axes.titlesize': 20,
                      'text.fontsize': 20,
                      'xtick.labelsize': 16,
                      'ytick.labelsize': 16,
                      'legend.pad': 0.2,     # empty space around the legend box
                      'legend.fontsize': 14,
                       'lines.markersize': 0,
                       'lines.linewidth': 3,
                      'font.size': 12,
                      'path.simplify': False}
        #              'figure.subplot.left':.10,
        #              'figure.subplot.bottom':.13,
        #              'figure.subplot.right':.94,
        #              'figure.subplot.top':.94}
        #              'figure.figsize': get_fig_size(800)}

        pylab.rcParams.update(plot_params)
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print "plotting ...."
        cax = ax.pcolormesh(d)#, edgecolor='k', linewidths='1')
        pylab.ylim(0, d.shape[0])
        pylab.xlim(0, d.shape[1])
        cbar = pylab.colorbar(cax)
        ax.set_xlabel('OR')
        ax.set_ylabel('Pattern number')
        cbar.set_label('Affinity')
          

        plot_fn = self.params['activation_matrix_fig']
        print "saving ....", plot_fn
        pylab.savefig(plot_fn, dpi=200)

        plot_fn = self.params['activation_matrix_fig'].rsplit('.png')[0] + '.pdf'
        print "saving ....", plot_fn
        pylab.savefig(plot_fn, dpi=200)


if __name__ == '__main__':

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        print 'Plotting default params'
        param_tool = simulation_parameters.parameter_storage()


    params = param_tool.params
    OAP = OdorantAffininityPlotter(params)
    OAP.plot_affinity_matrix()

pylab.show()
