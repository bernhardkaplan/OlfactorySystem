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
import MergeSpikefiles 
import simulation_parameters

class Plotter(object):
    def __init__(self, params, pn_max, cell_type):

        self.params = params
        Merger = MergeSpikefiles.MergeSpikefiles(params)

        d = np.zeros((pn_max, self.params['n_%s' % cell_type]))
        for pn in xrange(pn_max):
            print "Celltype: %s pattern_nr: %d" % (cell_type, pn)
            Merger.merge_nspike_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pn)
            Merger.merge_spiketimes_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pn)


            fn = self.params['%s_spikes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
            print "Loading data ", fn
            data = np.loadtxt(fn)
            idx = np.array(data[:, 0], dtype=np.int) - self.params['%s_offset' % cell_type]
#            print 'debug', idx.shape, d.shape, data.shape
            d[pn, idx] = data[:, 1]

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel('Cells')
        ax.set_ylabel('Pattern number')
        print "plotting ...."
        cax = ax.pcolormesh(d)#, edgecolor='k', linewidths='1')
        ax.set_xlim((0, d.shape[1]))
        ax.set_ylim((0, d.shape[0]))

        cbar = pylab.colorbar(cax)
        cbar.set_label('Number of spikes')
#            print "saving ...."
        pylab.show()

if __name__ == '__main__':

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params

    try:
        cell_type = sys.argv[2]
    except:
        print 'Missing cell_type argument'
        print 'Usage: python plot_activity_as_colormap.py FOLDER_NAME CELL_TYPE [PN_MAX]'
        exit(1)



    try:
        pn_max = int(sys.argv[3])
    except:
        print 'Plotting all patterns'
        pn_max = params['n_patterns']

#    pn_max = 4
    P = Plotter(params, pn_max, cell_type)


