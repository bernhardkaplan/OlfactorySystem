import os
import sys
import time
import numpy
import simulation_parameters # defines simulation parameters
import MergeSpikefiles
import MDSVQ
import AnalyseObOutput
# classes for setting up connectivity and the individual cell parameters
t1 = time.time()


if __name__ == '__main__':
    # ------------ I N I T -----------------------------
    # The simulation_parameters module defines a class for simulation parameter storage
    param_tool = simulation_parameters.parameter_storage()
    # params is the dictionary with all parameters
    params = param_tool.params
    param_tool.write_parameters_to_file(params["info_file"]) # human readable
    param_tool.write_parameters_to_file() # 


    # ------------ MDS + VQ of OB output ---------------
    ObAnalyser = AnalyseObOutput.AnalyseObOutput(params)
    ObAnalyser.get_output_file()
    ObAnalyser.get_output_activity()
    ObAnalyser.rescale_activity()


    ob_activity_fn = params['mit_response_normalized']

    mdsvq = MDSVQ.MDSVQ(params)
    activity_fn = params['mit_mds_input_fn'] # choose one of the different MIT - output files (with different normalizations)

    # 1) calculate mutual information (MI) between MIT cells and their distances
    # 2) run MDS in this MI-space (OB pattern reponse or mutual information (MI) space) and save mitral cell coordinates
    mds_output_fn = params['mds_ob_oc_output_fn']
    cell_type = 'mit'
    mdsvq.mds(activity_fn, mds_output_fn, thresh=1e-5, cell_type=cell_type)

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
    mdsvq.create_mitral_response_space(vq_output_fn, activity_fn, binary_oc_activation_fn, remove_silent_cells_fn=params['silent_mit_fn'], mit_mc_kmeans_trial=mit_mc_kmeans_trial) 

