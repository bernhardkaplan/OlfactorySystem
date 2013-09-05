import os 
import sys
import time
import MDSVQ
import BCPNN
import numpy as np
import AnalyseObOutput
import CreateObOcConnections

import CreatePyrReadoutParameters

def create_connectivity_after_mds_vq(params):
    PyrReadoutParamClass = CreatePyrReadoutParameters.CreatePyrReadoutParameters(params)
    PyrReadoutParamClass.write_pyr_parameters()
    PyrReadoutParamClass.write_readout_parameters()

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    print "Creating OB -> OC connections ..."
    ConnectionDrawer = CreateObOcConnections.CreateObOcConnections(param_tool, comm, rank, debug=1)

    if (rank == 0):
        ConnectionDrawer.get_mit_rsnp_connections_from_ALa(random_conn=params['ob_oc_random_conns'])
        ConnectionDrawer.get_mit_pyr_connections_from_ALa(random_conn=params['ob_oc_random_conns'])
        ConnectionDrawer.get_oc_oc_connections_from_ALa(random_conn=params['oc_oc_random_conns'])
        ConnectionDrawer.get_pyr_readout_connections_from_ALa()
    comm.Barrier()


    print "Creating hyper and minicolumns ..."
    ConnCreator = CreateHyperAndMinicolumns.CreateHyperAndMinicolumns(params)#, offset=0)
    ConnCreator.create_connections()
    comm.Barrier()

    #os.system("python select_gids_to_record.py")


def bcpnn_ob_oc(params):

    bcpnn = BCPNN.BCPNN(1, params['n_mit'], params['n_hc'], \
            params['n_mc'], params['n_patterns'], params) # 1 src HC, n_mit MC --> n_hc, n_mc

    binary_oc_activation_fn = params['binary_oc_activation_fn']
    w_ij_mit_hc = params['vq_ob_oc_output_fn']
    print "debug loading mask from:", w_ij_mit_hc

    activity_fn = params['oc_abstract_activity_fn']
    weights_fn = params['ob_oc_abstract_weights_fn']
    bias_fn = params['ob_oc_abstract_bias_fn']

    ob_activity_fn = params['mit_response_normalized']
    bcpnn.load_input_activity(ob_activity_fn)
    bcpnn.load_output_activity(binary_oc_activation_fn)
    bcpnn.load_mc_hc_mask(w_ij_mit_hc, silent_units_fn=params['silent_mit_fn'])
    bcpnn.initialize()

    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
    #    bcpnn.train_network()
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    #bcpnn.train_network(activity_fn, weights_fn, bias_fn)

    bcpnn.silence_mit(params['silent_mit_fn'])
    print "debug write to ", weights_fn
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)
    del bcpnn


def bcpnn_oc_oc(params):
    """
    train with the recurrent connections with OC output activity after learning OB -> OC
    """
    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, n_hc, n_mc, n_patterns, params)

    # train with binary oc activation derived from WTA after 2nd VQ
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    bcpnn.load_input_activity(binary_oc_activation_fn)
    bcpnn.load_output_activity(binary_oc_activation_fn)
    bcpnn.initialize()

    activity_fn = params['oc_rec_abstract_activity_fn'] # for output 
    weights_fn = params['oc_oc_abstract_weights_fn']
    bias_fn = params['oc_oc_abstract_bias_fn']
    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
    #    bcpnn.train_network()
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)
    del bcpnn


def bcpnn_oc_readout(params):

    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    n_readout = params['n_readout']
    if params['patterns_with_multiple_conc']:
        readout_activation = np.zeros((n_patterns, n_readout))
        for pn in xrange(n_patterns):
            active_readout_idx = pn / params['n_conc_per_pattern']
            readout_activation[pn, active_readout_idx] = 1
    else:
        readout_activation = np.eye(n_patterns)

    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    oc_activity_fn = params['oc_abstract_activity_fn']
    print 'Loading as input activity', oc_activity_fn
    bcpnn.load_input_activity(oc_activity_fn)
    bcpnn.load_output_activity(readout_activation)
    bcpnn.initialize()

    activity_fn = params['readout_abstract_activity_fn']
    weights_fn = params['oc_readout_abstract_weights_fn']
    bias_fn = params['oc_readout_abstract_bias_fn']

    #n_steps = 1
    n_steps = params['n_bcpnn_steps']
    for i in xrange(n_steps):
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    #bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)

    # testing
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    test_input = oc_activity_fn
    test_output = params['readout_abstract_activity_fn'].rsplit('.dat')[0] + '_test.dat'
    bcpnn.testing(test_input, weights_fn, bias_fn, test_output)
    del bcpnn



def mds_vq_ob_output(params):
    mdsvq = MDSVQ.MDSVQ(params)

    # 1) calculate mutual information (MI) between MIT cells and their distances
    # 2) run MDS in this MI-space (OB pattern reponse or mutual information (MI) space) and save mitral cell coordinates
    mds_output_fn = params['mds_ob_oc_output_fn']
    activity_fn = params['mit_mds_input_fn']
    mdsvq.mds(activity_fn, mds_output_fn, thresh=1e-5, cell_type='mit')

    # 3) run a VQ in the MI space, with n_clusters = n_hc
    # i.e. assign one or more hypercolumns to each mitral cell
    vq_output_fn = params['vq_ob_oc_output_fn']
    mdsvq.vq(mds_output_fn, vq_output_fn, params['n_hc'], overlap=params['vq_ob_oc_overlap'], remove_silent_cells_fn=params['silent_mit_fn'])

    # 4) For each hypercolumn, create a new space spanned by the mitral cells projecting to the hc
    #   Each mitral cell represents one dimension and each pattern represents one vector or point in that space.
    #   The value of one component in such a vector is equal to the normalized activation of the respective mitral cell.
    #   The n_patterns vectors are clustered by VQ among the minicolumns in the hypercolumn.
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    mdsvq.create_mitral_response_space(vq_output_fn, activity_fn, binary_oc_activation_fn, remove_silent_cells_fn=params['silent_mit_fn'])



if __name__ == '__main__':

    if len(sys.argv) > 1:
        param_fn = sys.argv[1]
        if os.path.isdir(param_fn):
            param_fn += '/Parameters/simulation_parameters.json'
        import json
        f = file(param_fn, 'r')
        print 'Loading parameters from', param_fn
        params = json.load(f)

    else:
        import simulation_parameters
        param_tool = simulation_parameters.parameter_storage()
        params = param_tool.params


    t_init = time.time()
    ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
    ObAnalyser.get_normalized_ob_activity()
    ObAnalyser.rescale_activity()

    t1 = time.time()
    mds_vq_ob_output(params)
    t1 = time.time()
    t2 = time.time()
    print "MDS-VQ Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)

    print 'BCPNN OB - OC'
    t1 = time.time()
    bcpnn_ob_oc(params)
    t2 = time.time()
    print "BCPNN-OB-OC Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)

    print 'BCPNN OC - OC'
    t1 = time.time()
    bcpnn_oc_oc(params)
    t2 = time.time()
    print "BCPNN OC - OC Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)

    print 'BCPNN OC - Readout'
    t1 = time.time()
    bcpnn_oc_readout(params)
    t2 = time.time()
    print "BCPNN OC - Readout %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)


#    mit_mc_kmeans_trial = 0
"""
import CreateOrnParameters
import CreateMitParameters
import CreatePyrReadoutParameters
import AnalyseObOutput
import CmapPlotter
import Optimizer
import SetOfCurvesPlotter
import RasterPlotEpthOb
import CreateSpikeTrain
import MergeSpikefiles
import GetConnectionsQuickFix
import CreateHyperAndMinicolumns
t1 = time.time()

# ------------ I N I T -----------------------------
# The network_parameters module defines a class for simulation parameter storage
#param_tool = network_parameters_BGL.simulation_parameters_BGL()
param_tool = network_parameters.simulation_parameters()
# params is the dictionary with all parameters
params = param_tool.params
param_tool.write_parameters_to_file(params["info_file"])
param_tool.print_cell_gids()

# write the simulation parameters to a NEURON executable file
folder_name = params['folder_name']

param_file_hoc = "%s/simulation_params.hoc" % folder_name
param_tool.hoc_file = "%s/simulation_params.hoc" % folder_name
param_tool.hoc_export()

#exit(1)

# EPTH -> OB connections are not affected by the pattern

#print "Creating connections: orn -> mit"
#ConnectionClass = CreateObConnections.CreateObConnections(params)
#ConnectionClass.connect_orn_mit()
#ConnectionClass.connect_orn_pg()
#ConnectionClass.connect_pg_mit_serial()
#ConnectionClass.connect_pg_mit_reciprocal()
#ConnectionClass.connect_mt_gran_local()
#ConnectionClass.connect_mt_gran_global()

#test = params["test"]
#OrnParamClass = CreateOrnParameters.CreateOrnParameters(params)
#if test == 1:
#    ok = OrnParamClass.create_test_parameters()
#else:
#    if (params['noisy_affinity_matrix']):
#        ok = OrnParamClass.create_parameters_orthogonal_patterns_new_rnd(params['OR_affinity_noise'])
#    else:
#        ok = OrnParamClass.create_parameters_orthogonal_patterns_new_rnd()
#if (ok != 1):
#    print "Something went wrong in OrnParamClass.create_test_parameters()"
#    exit(1)

#exit(1)

#MitParamClass = CreateMitParameters.CreateMitParameters(params)
#ok = MitParamClass.create_parameters(test=test)
#if (ok != 1):
#    print "Something went wrong in MitParamClass.create_parameters(n_patterns=n_patterns"
#    exit(1)


#ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
#ObAnalyser.get_output_activity()
#ObAnalyser.rescale_activity()
#ObAnalyser.rescale_activity_cellwise()
#ObAnalyser.rescale_activity_glom_patterns()
#ObAnalyser.rescale_activity_patternwise()
#ObAnalyser.rescale_activity_patternwise_then_cellwise()
#ObAnalyser.rescale_activity_cellwise_then_patternwise()
#ObAnalyser.get_output_file()

#os.system("python bcpnn_all.py")

PyrReadoutParamClass = CreatePyrReadoutParameters.CreatePyrReadoutParameters(params)
PyrReadoutParamClass.write_pyr_parameters()
PyrReadoutParamClass.write_readout_parameters()

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print "Creating OB -> OC connections ..."
GetConn = GetConnectionsQuickFix.GetConnectionsQuickFix(param_tool, comm, rank, debug=1)

if (rank == 0):
    GetConn.get_mit_rsnp_connections_from_ALa(random_conn=params['ob_oc_random_conns'])
    GetConn.get_mit_pyr_connections_from_ALa(random_conn=params['ob_oc_random_conns'])
    GetConn.get_oc_oc_connections_from_ALa(random_conn=params['oc_oc_random_conns'])
    GetConn.get_pyr_readout_connections_from_ALa()
comm.Barrier()


print "Creating hyper and minicolumns ..."
ConnCreator = CreateHyperAndMinicolumns.CreateHyperAndMinicolumns(params)#, offset=0)
ConnCreator.create_connections()
comm.Barrier()

#os.system("python select_gids_to_record.py")

print "Folder name:", params['folder_name']
t2 = time.time()
print "Time: %.1f sec %.1f min" % (t2 - t1, (t2-t1)/60.)
"""
