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
    def __init__(self, params):

        self.params = params
        self.Merger = MergeSpikefiles.MergeSpikefiles(params)
        self.fontsize = 18

    def plot(self, pn_max, cell_type, mc_clustering=False):
        self.cell_type = cell_type
        self.pn_max = pn_max
        d = np.zeros((pn_max, self.params['n_%s' % cell_type]))
        for pn in xrange(pn_max):

            print "Celltype: %s pattern_nr: %d" % (cell_type, pn)
            fn = self.params['%s_spikes_merged_fn_base' % (cell_type)] + str(pn) + '.dat'
            if not os.path.exists(fn):
                print 'DEBUG Merging ...\n\n'
                self.Merger.merge_nspike_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pn)
                self.Merger.merge_spiketimes_files(self.params['%s_spike_fn_base' % (cell_type)], self.params['%s_spikes_merged_fn_base' % (cell_type)], pn)

            print "Loading data ", fn
            data = np.loadtxt(fn)
            idx = np.array(data[:, 0], dtype=np.int) - self.params['%s_offset' % cell_type]
#            print 'debug', idx.shape, d.shape, data.shape
            d[pn, idx] = data[:, 1]

        if mc_clustering:
            d = self.cluster_nspikes_by_mc(d)
            xlabel = 'Minicolumn'
            clabel = 'Average number of spikes'
            fig_fn = self.params['figure_folder'] + '/' + '%s_activity_%dpatterns_mcclustered.png' % (self.cell_type, self.pn_max)
        else:
            xlabel = 'Cells'
            clabel = 'Number of spikes'
            fig_fn = self.params['figure_folder'] + '/' + '%s_activity_%dpatterns.png' % (self.cell_type, self.pn_max)

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.set_title('Activity for %s cells over %d patterns' % (cell_type.capitalize(), pn_max), fontsize=self.fontsize)
        ax.set_xlabel(xlabel, fontsize=self.fontsize)
        ax.set_ylabel('Pattern number', fontsize=self.fontsize)
        print "plotting ...."
        cax = ax.pcolormesh(d)#, cmap='binary')
        ax.set_xlim((0, d.shape[1]))
        ax.set_ylim((0, d.shape[0]))

        cbar = pylab.colorbar(cax)
        cbar.set_label(clabel, fontsize=self.fontsize)
        print 'Saving to:', fig_fn
        pylab.savefig(fig_fn)
        pylab.show()


    def cluster_nspikes_by_mc(self, all_nspikes):

        n_col = self.params['n_hc'] * self.params['n_mc']
        d_out = np.zeros((all_nspikes.shape[0], n_col))
        n_cells_per_mc = self.params['n_%s_per_mc' % self.cell_type]
        for pn in xrange(self.pn_max):
            for mc in xrange(n_col):
                idx0 = mc * n_cells_per_mc
                idx1 = (mc + 1) * n_cells_per_mc
                d_out[pn, mc] = all_nspikes[pn, idx0:idx1].sum() / n_cells_per_mc

        return d_out


if __name__ == '__main__':

    try:
        folder = sys.argv[1]
        params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
        param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
    except:
        param_tool = simulation_parameters.parameter_storage()

    params = param_tool.params
    print 'debug n_cells', params['n_cells']

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
    P = Plotter(params)

    if (cell_type == 'pyr' or cell_type == 'rsnp'):
        ok = raw_input('Cluster nspikes by minicolumn?\n')
        if ok.capitalize() == 'Y' or ok == '':
            P.plot(pn_max, cell_type, mc_clustering=True)
        else:
            P.plot(pn_max, cell_type)
    else:
        P.plot(pn_max, cell_type)
