"""

    To run this script make sure that in simulation_parameters.py:
    n_or = 32
    n_gor = number_of_lists ( = number of curves)
    rel_orn_mit = 1

"""

import simulation_parameters
import os
import time
import MergeSpikefiles
import SetOfCurvesPlotter
import CreateOrnParameters
import numpy as np
import fit_orn_params

param_tool = simulation_parameters.parameter_storage()
params = param_tool.params
param_tool.hoc_export()

# class to write parameters into file
OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
OrnParamClass.param_list = ["gna", "gk", "gkcag", "gcal", "gleak_orn", "tau_cadec"]

# define parameters:
"""
each column gets a own list of parameters
this list is passed to OrnParamClass via the member current_param_values
e.g.
col = 0
OrnParamClass.current_param_values[col][2] 
corresponds to the "gkcag" values for all ORNs in this column
"""

OrnParamClass.set_oor_value_for_conc_sweep()
OrnParamClass.current_param_values = np.zeros((params["n_orn_x"], len(OrnParamClass.param_list)))#.tolist()

#min_params =    [0.5, 5e-2, 1.3e-3,     1.3e-4,     7.5e-5, 1000]
#list2 =         [0.5, 5e-2, 1.35e-3,    1.35e-4,    7.4e-5, 1000]
#list3 =         [0.5, 5e-2, 2.0e-3,     2.0e-4,     7.3e-5, 1000]
#list4 =         [0.5, 5e-2, 3.25e-3,    3.25e-4,    7.2e-5, 1000]
#list5 =         [0.5, 5e-2, 4.3e-3,     4.3e-4,     7.1e-5, 1000]
#list6 =         [0.5, 5e-2, 5.5e-3,     5.5e-4,     7.0e-5, 1000]
#list7 =         [0.5, 5e-2, 6.5e-3,     6.5e-4,     7.0e-5, 1000]
#list8 =         [0.5, 5e-2, 7.5e-3,     7.5e-4,     6.0e-5, 1000]
#list9 =         [0.5, 5e-2, 8.5e-3,     8.5e-4,     5.0e-5, 1000]
#list10=         [0.5, 5e-2, 9.5e-3,     9.5e-4,     4.0e-5, 1000]
#list11 =        [0.5, 5e-2, 9e-3,       9e-4,       3.0e-5, 1000]
#list12 =        [0.5, 5e-2, 1e-2,       1e-3,       2.0e-5, 1000]
#list13 =        [0.5, 5e-2, 1.1e-2,     1.1e-3,     1.0e-5, 1000]
#list14 =        [0.5, 5e-2, 1.2e-2,     1.2e-3,     9.0e-4, 1000]
#list15 =        [0.5, 5e-2, 1.3e-2,     1.3e-3,     8.0e-4, 1000]
#max_params =    [0.5, 5e-2, 1.4e-2,     1.4e-3,     7.0e-4, 1000]

#for param in xrange(len(OrnParamClass.param_list)):
#    OrnParamClass.current_param_values[0][param] = min_params[param]
#    OrnParamClass.current_param_values[1][param] = list2[param]
#    OrnParamClass.current_param_values[2][param] = list3[param]
#    OrnParamClass.current_param_values[3][param] = list4[param]
#    OrnParamClass.current_param_values[4][param] = list5[param]
#    OrnParamClass.current_param_values[5][param] = list6[param]
#    OrnParamClass.current_param_values[6][param] = list7[param]
#    OrnParamClass.current_param_values[7][param] = list8[param]
#    OrnParamClass.current_param_values[8][param] = list9[param]
#    OrnParamClass.current_param_values[9][param] = list10[param]
#    OrnParamClass.current_param_values[10][param] = list11[param]
#    OrnParamClass.current_param_values[11][param] = list12[param]
#    OrnParamClass.current_param_values[12][param] = list13[param]
#    OrnParamClass.current_param_values[13][param] = list14[param]
#    OrnParamClass.current_param_values[14][param] = list15[param]
#    OrnParamClass.current_param_values[15][param] = max_params[param]


#list1 =         [0.5, 5e-2, 1.3e-3, 1.3e-4,     7.75e-5, 1000]
#list2 =         [0.5, 5e-2, 1.35e-3, 1.35e-4,     7.7e-5, 1000]
#list3 =         [0.5, 5e-2, 2.0e-3, 2.0e-4,     7.675e-5, 1000]
#list4 =         [0.5, 5e-2, 3.25e-3, 3.250e-4,     7.65e-5, 1000]
#list5 =         [0.5, 5e-2, 4.3e-3, 4.3e-4,     7.6e-5, 1000]
#list6 =         [0.5, 5e-2, 5.5e-3, 5.5e-4,     7.5e-5, 1000]

list1 =         [0.5, 5e-2, 1.3e-3, 1.3e-4,     7.5e-5, 1000]
list2 =         [0.5, 5e-2, 1.35e-3, 1.35e-4,   7.4e-5, 1000]
list3 =         [0.5, 5e-2, 2.0e-3, 2.0e-4,     7.3e-5, 1000]
list4 =         [0.5, 5e-2, 3.25e-3, 3.250e-4,  7.2e-5, 1000]
list5 =         [0.5, 5e-2, 4.3e-3, 4.3e-4,     7.1e-5, 1000]
list6 =         [0.5, 5e-2, 5.5e-3, 5.5e-4,     7.0e-5, 1000]
list7 =         [0.5, 5e-2, 6.5e-3, 6.5e-4,     7.0e-5, 1000]
list8 =         [0.5, 5e-2, 7.5e-3, 7.5e-4,     6.0e-5, 1000]
list9 =         [0.5, 5e-2, 8.5e-3, 8.5e-4,     5.0e-5, 1000]
list10=         [0.5, 5e-2, 9.5e-3, 9.5e-4,     4.0e-5, 1000]
for param in xrange(len(OrnParamClass.param_list)):
    OrnParamClass.current_param_values[0][param] = list1[param]
    OrnParamClass.current_param_values[1][param] = list2[param]
    OrnParamClass.current_param_values[2][param] = list3[param]
    OrnParamClass.current_param_values[3][param] = list4[param]
    OrnParamClass.current_param_values[4][param] = list5[param]
    OrnParamClass.current_param_values[5][param] = list6[param]
    OrnParamClass.current_param_values[6][param] = list7[param]
    OrnParamClass.current_param_values[7][param] = list8[param]
    OrnParamClass.current_param_values[8][param] = list9[param]
    OrnParamClass.current_param_values[9][param] = list10[param]


#min_params =    [0.5, 5e-2, 1e-4, 1e-5, 2.1e-5, 1000]
#list2 =         [0.5, 5e-2, 1e-4, 1e-5, 2.1e-5, 1000]
#list3 =         [0.5, 5e-2, 1e-4, 1e-5, 2.55e-5, 1000]
#list4 =         [0.5, 5e-2, 1e-3, 1e-4, 4e-5, 1000]
#list5 =         [0.5, 5e-2, 2e-3, 2e-4, 6e-5, 1000]
#list6 =         [0.5, 5e-2, 3.15e-3, 3.15e-4, 8e-5, 1000]
#list7 =         [0.5, 5e-2, 4.25e-3, 4.25e-4, 1.1e-4, 1000]
#list8 =         [0.5, 5e-2, 5.5e-3, 5.5e-4, 1.3e-4, 1000]
#list9 =         [0.5, 5e-2, 7e-3, 7e-4, 1.3e-4, 1000]
#list10 =         [0.5, 5e-2, 8e-3, 8e-4, 1.3e-4, 1000]
#list11 =         [0.5, 5e-2, 9e-3,  9e-4, 1.3e-4, 1000]
#list12 =         [0.5, 5e-2, 1e-2,  1e-3, 1.3e-4, 1000]
#list13 =         [0.5, 5e-2, 1.1e-2, 1.1e-3, 1e-4, 1000]
#list14 =         [0.5, 5e-2, 1.2e-2, 1.2e-3, 1e-4, 1000]
#list15 =         [0.5, 5e-2, 1.3e-2, 1.3e-3, 5e-5, 1000]
#max_params =    [0.5, 5e-2, 1.4e-2, 1.4e-3, 1e-5, 1000]
#for param in xrange(len(OrnParamClass.param_list)):
#    OrnParamClass.current_param_values[0][param] = min_params[param]
#    OrnParamClass.current_param_values[1][param] = list2[param]
#    OrnParamClass.current_param_values[2][param] = list3[param]
#    OrnParamClass.current_param_values[3][param] = list4[param]
#    OrnParamClass.current_param_values[4][param] = list5[param]
#    OrnParamClass.current_param_values[5][param] = list6[param]
#    OrnParamClass.current_param_values[6][param] = list7[param]
#    OrnParamClass.current_param_values[7][param] = list8[param]
#    OrnParamClass.current_param_values[8][param] = list9[param]
#    OrnParamClass.current_param_values[9][param] = list10[param]
#    OrnParamClass.current_param_values[10][param] = list11[param]
#    OrnParamClass.current_param_values[11][param] = list12[param]
#    OrnParamClass.current_param_values[12][param] = list13[param]
#    OrnParamClass.current_param_values[13][param] = list14[param]
#    OrnParamClass.current_param_values[14][param] = list15[param]
#    OrnParamClass.current_param_values[15][param] = max_params[param]


        
# now overwrite gleak values according to the corresponding function
#OrnParamClass.set_gleak_values()
OrnParamClass.write_current_param_values_to_file()

os.chdir('neuron_files') # this is important to avoide problems with tabchannel files and the functions defined therein

pn = 0
neuron_command = "mpirun -np %d $(which nrniv) -mpi -nobanner -nogui \
        -c \"x=%d\" -c \"strdef param_file\" -c \"sprint(param_file, \\\"%s\\\")\" start_file_epth_response_curve.hoc > delme%d" \
        % (params['n_proc'], pn, params['hoc_file'], pn)

os.system("rm %s/*" % (params["spiketimes_folder"]))
os.system("rm %s/*" % (params["volt_folder"]))
os.system(neuron_command)



Merger = MergeSpikefiles.MergeSpikefiles(params)
Merger.merge_epth_spiketimes_file(pattern=pn)
Merger.merge_epth_nspike_files(pattern=pn)

SOCP = SetOfCurvesPlotter.SetOfCurvesPlotter(params)

sim_cnt = 2
log_file = open("params_handtuning_%d.txt" % sim_cnt, 'w')
for i in xrange(params["n_gor"]):
    line = "%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e\n" % ( \
            OrnParamClass.gor_values[i],\
            OrnParamClass.current_param_values[i][0],\
            OrnParamClass.current_param_values[i][1],\
            OrnParamClass.current_param_values[i][2],\
            OrnParamClass.current_param_values[i][3],\
            OrnParamClass.current_param_values[i][4],\
            OrnParamClass.current_param_values[i][5])
    log_file.write(line)
    log_file.flush()

output_fn = params['figure_folder'] + '/hand_tuned_%d.png' % sim_cnt
print 'Saving figure to:', output_fn
SOCP.plot_set_of_curves()
#SOCP.plot_set_of_curves(output_fn=output_fn)
