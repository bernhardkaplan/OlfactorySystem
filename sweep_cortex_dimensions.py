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
import numpy as np
import simulation_parameters # defines simulation parameters
import prepare_epth_ob_prelearning
import MergeSpikefiles
import MDSVQ
import BCPNN
import AnalyseObOutput
import CreatePyrReadoutParameters
import CreateHyperAndMinicolumns
import TransformAbstractToDetailedConnectivity
import CreateObOcConnections as ObOc
# classes for setting up connectivity and the individual cell parameters
t_init = time.time()


def check_readout_test(fn):
    d = np.loadtxt(fn)
    n_correct = 0
    for row in xrange(d.shape[0]):
        winner = d[row, :].argmax()
        if row == winner:
            n_correct += 1
    return n_correct




if __name__ == '__main__':


    # all setups should use the same MI-coordinates for the MT cells
    mds_output_fn = 'Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/Other/mds_ob_oc_output.dat'
    mit_response_normalized_fn = 'Cluster_SparserObPatterns_nGlom40_nHC9_nMC9_vqOvrlp8_ORnoise0.0_OrAffNorm0_postL_np50_1_OcOnly/NumberOfSpikes/mit_response_normalized_np50'
    n_hc_mc = [(6, 60), (8, 45), (12, 30), (18,20), (20, 18), (30, 12), (45, 8), (60, 6)]

    log_file = open('ctx_params_sweep.txt', 'w')
#    n_hc_mc = [(6, 60)]
    for i_ in xrange(len(n_hc_mc)):
        n_hc, n_mc = n_hc_mc[i_][0], n_hc_mc[i_][1]
#        for vq_overlap in [0]:
        for vq_overlap in [0, 1, 2, 3, 4, 5, n_hc - 1]:

            print 'n_hc, n_mc', n_hc, n_mc

            param_tool = simulation_parameters.parameter_storage()
            params = param_tool.params
            params['n_hc'] = n_hc
            params['n_mc'] = n_mc
            params['vq_ob_oc_overlap'] = vq_overlap
            param_tool.set_filenames()

            cmd = 'cp %s %s' % (mit_response_normalized_fn, params['nspikes_folder'])
            print 'Copying', cmd
            os.system(cmd)

            ObOc.mds_vq_ob_output(params, mds_output_fn)
            ObOc.bcpnn_ob_oc(params)
            ObOc.bcpnn_oc_oc(params)
            ObOc.bcpnn_oc_readout(params)
            
            test_output = params['readout_abstract_activity_fn'].rsplit('.dat')[0] + '_test.dat'
#            wta_output_fn = output_fn.rsplit('.dat')[0] + '_wta.dat'
            n_correct = check_readout_test(test_output)
            log_txt = '%d\t%d\t%d\t%d\t%.2f\n' % (n_hc, n_mc, vq_overlap, n_correct, n_correct / float(params['n_patterns']))
            print 'log_txt:', log_txt
            log_file.write(log_txt)
            log_file.flush()

    log_file.close()


