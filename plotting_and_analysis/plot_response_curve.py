import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
print 'cmd_subfolder', cmd_subfolder
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import simulation_parameters
import MergeSpikefiles
import SetOfCurvesPlotter

info_txt = \
"""
Usage:
    python plot_response_curve.py [FOLDER] [CELLTYPE]
"""
assert (len(sys.argv) > 2), 'ERROR: folder and cell_type not given\n' + info_txt
folder = sys.argv[1]
cell_type = sys.argv[2]

params_fn = os.path.abspath(folder) + '/Parameters/simulation_parameters.json'
param_tool = simulation_parameters.parameter_storage(params_fn=params_fn)
params = param_tool.params
print 'debug', params['%s_spikes_merged_fn_base' % cell_type]

pn = 0
Merger = MergeSpikefiles.MergeSpikefiles(params)
if cell_type.lower() == 'mit':
    Merger.merge_ob_spiketimes_file(pattern=pn)
    Merger.merge_ob_nspike_files(pattern=pn)
elif cell_type.lower() == 'orn':
    Merger.merge_epth_spiketimes_file(pattern=pn)
    Merger.merge_epth_nspike_files(pattern=pn)
else:
    print 'ERROR: Invalid cell type given %s\n' % cell_type


sim_cnt = 0
SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)
output_fn = params['figure_folder'] + '/ob_response_curve_%d.png' % sim_cnt
SOCP.plot_set_of_curves(output_fn, cell_type='mit')
print 'Opening with ristretto: %s' % (output_fn)
os.system('ristretto %s' % output_fn)

# ------- Merge spike files
