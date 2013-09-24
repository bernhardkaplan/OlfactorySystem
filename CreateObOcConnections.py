info_txt = \
"""
    Before running this script, you should have:
    - the simulation results of a 'pre-learning' run 

    - chosen a new folder name and run this script 
    only creating the new folder structure and parameter files by running (see below)
        param_tool.write_parameters_to_file() # 
    - then, copy the nspike files of mitral cells into the new NumberOfSpikes folder:
     cp PRE_LEARNING/NumberOfSpikes/mit*merged* POST_LEARNING/NumberOfSpikes/
"""

import os
import sys
import time
import numpy
import simulation_parameters # defines simulation parameters
import prepare_epth_ob_prelearning
import MergeSpikefiles
import MDSVQ
import BCPNN
import AnalyseObOutput
import CreatePyrReadoutParameters
import CreateHyperAndMinicolumns
import TransformAbstractToDetailedConnectivity
# classes for setting up connectivity and the individual cell parameters
t_init = time.time()

def mds_vq_ob_output(params):
    ob_activity_fn = params['mit_response_normalized']

    mdsvq = MDSVQ.MDSVQ(params)
    activity_fn = params['mit_mds_input_fn'] # choose one of the different MIT - output files (with different normalizations)

    # 1) calculate mutual information (MI) between MIT cells and their distances
    # 2) run MDS in this MI-space (OB pattern reponse or mutual information (MI) space) and save mitral cell coordinates
    mds_output_fn = params['mds_ob_oc_output_fn']
    cell_type = 'mit'
    mdsvq.mds(activity_fn, mds_output_fn, thresh=1e-6, cell_type=cell_type)

    t_2 = time.time()
    t_diff = t_2 - t_init
    print "time in sec for mds:", t_diff

    # 3) run a VQ in the MI space, with n_clusters = n_hc
    # i.e. assign one or more hypercolumns to each mitral cell
    vq_output_fn = params['vq_ob_oc_output_fn']
    mdsvq.vq(mds_output_fn, vq_output_fn, params['n_hc'], overlap=params['vq_ob_oc_overlap'], remove_silent_cells_fn=params['silent_mit_fn'])
    t_3 = time.time()
    t_diff = t_3 - t_2
    print "time in sec for vq:", t_diff 

    # 4) For each hypercolumn, create a new space spanned by the mitral cells projecting to the hc
    #   Each mitral cell represents one dimension and each pattern represents one vector or point in that space.
    #   The value of one component in such a vector is equal to the normalized activation of the respective mitral cell.
    #   The n_patterns vectors are clustered by VQ among the minicolumns in the hypercolumn.
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    mit_mc_kmeans_trial = 0
    mdsvq.create_mitral_response_space(vq_output_fn, activity_fn, binary_oc_activation_fn, remove_silent_cells_fn=params['silent_mit_fn'], mit_mc_kmeans_trial=mit_mc_kmeans_trial) 
    return (ob_activity_fn, binary_oc_activation_fn, vq_output_fn)


def bcpnn_ob_oc(params):

    n_mit = params['n_mit']

    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    bcpnn = BCPNN.BCPNN(1, n_mit, n_hc, n_mc, n_patterns, params) # 1 src HC, n_mit MC --> n_hc, n_mc

    ob_activity_fn = params['mit_response_normalized']
    activity_fn = params['oc_abstract_activity_fn']
    
    bias_fn = params['ob_oc_abstract_bias_fn']
    binary_oc_activation_fn = params['binary_oc_activation_fn']
    w_ij_mit_hc = params['vq_ob_oc_output_fn']
    weights_fn = params['ob_oc_abstract_weights_fn']
    print "BCPNN OB -> OC loading files:", ob_activity_fn, '\n', activity_fn, '\n', weights_fn, '\n', bias_fn, '\n', binary_oc_activation_fn, '\n', w_ij_mit_hc


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
#    print 'bcpnn_ob_oc: input = %s\tweights = %s\tbias = %s\ttest_output = %s' % (test_input, weights_fn, bias_fn, test_output)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)
    del bcpnn

def bcpnn_oc_oc(params):

    n_patterns = params['n_patterns']
    n_hc = params['n_hc']
    n_mc = params['n_mc']
    n_readout = params['n_readout']
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, n_hc, n_mc, n_patterns, params)

    # train with the recurrent connections with OC output activity after learning OB -> OC
    #oc_activity_fn = params['oc_abstract_activity_fn']
    #bcpnn.load_input_activity(oc_activity_fn)
    #bcpnn.load_output_activity(oc_activity_fn)
    #bcpnn.initialize()

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
#    if params['multiple_concentrations_per_pattern']:
#        n_readout = self.params['SOMETHING_TODO']
#        readout_activation = numpy.zeros((n_patterns, n_readout))
#        for pn in xrange(n_patterns):
#            active_readout_idx = pn / params['n_concentrations']
#            readout_activation[pn, active_readout_idx] = 1
#    else:
    n_readout = params['n_patterns']
    readout_activation = numpy.eye(n_patterns)
    print 'readout_activation', readout_activation

    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    #oc_activity_fn = params['oc_abstract_activity_fn'].rsplit('.dat')[0] + '0.dat'
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
    #    bcpnn.train_network()
        bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    #bcpnn.train_network(activity_fn, weights_fn, bias_fn)
    bcpnn.write_to_files(activity_fn, weights_fn, bias_fn)


    #if params['multiple_concentrations_per_pattern']:
    #    n_patterns = 10
    # testing 
    #del bcpnn
    bcpnn = BCPNN.BCPNN(n_hc, n_mc, 1, n_readout, n_patterns, params)
    test_input = oc_activity_fn
    test_output = params['readout_abstract_activity_fn'].rsplit('.dat')[0] + '_test.dat'
    print 'BCPNN.testing(input=%s, weights=%s, bias=%s, test_output=%s' % (test_input, weights_fn, bias_fn, test_output)
    bcpnn.testing(test_input, weights_fn, bias_fn, test_output)



def create_pyr_parameters(params):

    PyrReadoutParamClass = CreatePyrReadoutParameters.CreatePyrReadoutParameters(params)
    PyrReadoutParamClass.write_pyr_parameters()
    PyrReadoutParamClass.write_readout_parameters()


def create_connections(params):
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    print "Creating OB -> OC connections ..."
    GetConn = TransformAbstractToDetailedConnectivity.GetConnectionsQuickFix(param_tool, comm, rank, debug=1)

    if (rank == 0):
        GetConn.get_mit_rsnp_connections(random_conn=params['ob_oc_random_conns'])
        GetConn.get_mit_pyr_connections(random_conn=params['ob_oc_random_conns'])
        GetConn.get_oc_oc_connections(random_conn=params['oc_oc_random_conns'])
        GetConn.get_pyr_readout_connections()
    comm.Barrier()

    print "Creating hyper and minicolumns ..."
    ConnCreator = CreateHyperAndMinicolumns.CreateHyperAndMinicolumns(params)#, offset=0)
    ConnCreator.create_connections()
    comm.Barrier()

    print "Folder name:", params['folder_name']


if __name__ == '__main__':
    print info_txt
    # ------------ I N I T -----------------------------
    # The simulation_parameters module defines a class for simulation parameter storage
    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params
    param_tool.write_parameters_to_file(params["info_file"])
    param_tool.write_parameters_to_file() # 
    param_tool.hoc_export() # 


#    prepare_epth_ob_prelearning.prepare_epth_ob(params)

#     ------------ MDS + VQ of OB output ---------------
#    ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
#    ObAnalyser.get_output_file()
#    ObAnalyser.get_output_activity()
#    ObAnalyser.rescale_activity()
#    mds_vq_ob_output(params)

#    bcpnn_ob_oc(params)
#    bcpnn_oc_oc(params)
#    bcpnn_oc_readout(params)

#    create_pyr_parameters(params)
#    create_connections(params)

