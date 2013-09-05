import sys
import re
import numpy as np
import numpy.random as rnd
import random
import os
import pylab
import time
import matplotlib
import utils

class CreateObOcConnections(object):
    def __init__(self, param_tool, comm=None, rank=0, debug=0):
        '''
        param_tool: class storing the dictionary with all simulation parameters
        comm is the MPI communicator
        rank is the mpi rank of the process
        '''
        self.param_tool = param_tool
        self.params = param_tool.params
        self.comm = comm    # MPI.COMM_WORLD
        self.debug = debug
        if comm != None:
            self.n_proc = self.comm.Get_size()
        else:
            self.n_proc = 1
        self.my_rank = rank # process id


    def check_vq_and_conn_files(self):
        '''
        check the dimensionality of the vq file and the hebbconnection_ob_oc files
        each line in the vq file contains the 
        '''

        glom_hc_mapping = np.ones(n_mit_y) # to which tgt_hc in cortex does the src_hc (glom) connect to?
        glom_hc_mapping *= magic_number
        # glom_hc_mapping[src_hc] = tgt_hc
        for tgt_hc in xrange(num_tgt_hc): #tgt_hc is equivalent to a line in vq_file
            # remove whitespace from the vq lines
            src_hcs = re.sub('\s', '', vq_lines[tgt_hc]) #src_hcs = '2,7'
            if (src_hcs != ''): # if this tgt_hc gets input and is not an empty line in the vq_file ...
                # get the list of source hypercolumns
                src_hcs = src_hcs.rsplit(',') # e.g. src_hcs = [2, 7] 
                for src_hc in src_hcs:
                    glom_hc_mapping[int(src_hc)] = int(tgt_hc)
                    vq_data[tgt_hc].append(int(src_hc))


    def get_vq_result(self):
        """
        Returns the result of the VQ as a matrix in the following format:
            row = glomerulus index
            column = target Hypercolumn
        if conn_mat[glom, hc] == 1:
            glom projects to hc
        else:
            no connection between glom and hc
            i.e.
            conn_mat[glom, hc] = 0 
        """

        # get data from vq_file
        vq_file = file(self.params["vq_path"], "r")
        vq_data = []
        for line in vq_file:
            vq_data.append(line.rsplit(','))
        num_tgt_hc = len(vq_data)
        n_glom = self.params['n_mit_y']
        n_hc = self.params['n_hc']
        conn_mat = np.zeros((n_glom, n_hc))
        assert (len(vq_data) == n_hc), "Number of HC in network_parameters does not match the number of HC in the vq_file to plot"
        for hc in xrange(len(vq_data)):
            for i in xrange(len(vq_data[hc])):
                src_glom = vq_data[hc][i]
                conn_mat[src_glom, hc] = 1.

        return conn_mat


    def get_abstract_weight_matrix_ob_oc(self):
        """
        Returns the abstract weights obtained from the vq file and the hebbconnection file as two matrices:
            abstract_weights (with raw values, i.e. unscaled before taking the log)
            abstract_weights_log (after scaling the weights with params['ob_oc_weight_scaling'] and then taking the log)
        Additionally, these two matrices are stored in the two files named:
         self.params['connection_matrix_abstract_ob_oc_dat']
         self.params['connection_matrix_abstract_ob_oc_dat_unscaled']

        """
        vq_file = file(self.params["vq_path"], "r") # get data from vq_file
        vq_data = []
        for line in vq_file:
            vq_data.append(line.rsplit(','))
        num_tgt_hc = len(vq_data)
        # ------------- get abstract weights ------------------ 
        # read lines from the Conn_file which stores the weights for
        # tgt_mc vs src_mc (ids for the hcs is stored in vq_file)
        ob_oc_abstract_conns = file(self.params['mit_pyr_conn_path'], 'r')
        conn_file_lines = ob_oc_abstract_conns.readlines()
        mit_mc_conns = [] # a list of lists with conn_file_lines[tgt_mc] = src_mc

        tgt_mc = 0
        for line_cnt in xrange(len(conn_file_lines)): #line_cnt is a tgt_mc
            # remove all the whitespace
            conn_file_lines[line_cnt] = re.sub("\s", "", conn_file_lines[line_cnt])
            mit_mc_conns.append(conn_file_lines[line_cnt].rsplit(","))
            tgt_mc += 1
        n_mit = self.params["n_mit"]
        n_mit_x = self.params["n_mit_x"]
        n_pyr = self.params["n_pyr"]
        n_mc = self.params["n_mc"]
        n_hc = self.params["n_hc"]
        abstract_weights = np.zeros((n_mit, n_hc * n_mc)) # abstract weights
        abstract_weights_log = np.zeros((n_mit, n_hc * n_mc)) # abstract weights

        w_min = 0.0
        w_max = 0.0
        assert (len(mit_mc_conns) == n_hc * n_mc)
        for tgt_mc in xrange(n_hc * n_mc): # iterate through rows in hebbconnections_ob_oc.txt 
            for w_cnt in xrange(len(mit_mc_conns[tgt_mc])): # for all given weights
                # calculate src index
                tgt_hc = tgt_mc / n_mc # = line_in_vq_file
                element_in_vq_line = w_cnt / n_mit_x
                src_glom = int(vq_data[tgt_hc][element_in_vq_line])
                src_index = src_glom * n_mit_x + w_cnt % n_mit_x
                tgt_index = tgt_hc * n_mc + tgt_mc % n_mc
                w = float(mit_mc_conns[tgt_mc][w_cnt])

                if (w != 0.):
                    w_log = np.log(w * self.params['ob_oc_weight_scaling'])
                    if (w_log) > 0:
                        w_max = max(w_log, w_max)
                    else:
                        w_min = min(w_log, w_min)
                    abstract_weights_log[src_index, tgt_index] = w_log
                    
                abstract_weights[src_index, tgt_index] = w
    
        output_fn = self.params['connection_matrix_abstract_ob_oc_dat']
        np.savetxt(output_fn, abstract_weights_log)
        output_fn_unscaled = self.params['connection_matrix_abstract_ob_oc_dat_unscaled']
        np.savetxt(output_fn_unscaled, abstract_weights)
        return abstract_weights, abstract_weights_log



    def get_mit_pyr_connections_from_ALa(self, random_conn=False):
        '''
        This functions creates the OB->OC connections files from the matrix stored in self.params['ob_oc_abstract_weights_fn']
        '''
        print "Drawing MIT - PYR connections .... "
        abstract_weights = np.loadtxt(self.params['ob_oc_abstract_weights_fn'])
        if random_conn:
            rnd.seed(self.params['random_ob_oc_seed'])
            rnd.shuffle(abstract_weights)
            np.savetxt(self.params['ob_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)
        assert (abstract_weights[:, 0].size == self.params['n_mit'])
        assert (abstract_weights[0, :].size == self.params['n_hc'] * self.params['n_mc'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()
        w_mit_pyr_max = self.params['w_mit_pyr_max']
        w_mit_pyr_matrix = np.zeros((self.params['n_mit'], self.params['n_hc'] * self.params['n_mc']))

        output = ""
        line_cnt = 0
        for tgt_mc in xrange(self.params['n_hc'] * self.params['n_mc']): #column
            for mit in xrange(self.params['n_mit']): # row
                w_in = abstract_weights[mit, tgt_mc]
                if (w_in > 0):
#                    for tgt_pyr in xrange(int(round(self.params['n_tgt_pyr_per_mc']))):
                    tgt_pyrs = utils.get_rnd_targets(self.params['n_pyr_per_mc'], self.params['n_tgt_pyr_per_mc'])
                    for tgt_pyr in tgt_pyrs:
                        w_out = (w_in / w_max_abstract) * w_mit_pyr_max
                        w_noise = rnd.normal(w_out, w_out * self.params['w_mit_pyr_sigma_frac'])
                        if (w_noise > self.params['weight_threshold']):
                            src_gid = self.params['mit_offset'] + mit
                            tgt_gid = self.params['pyr_offset'] + tgt_pyr + tgt_mc * self.params['n_pyr_per_mc']
                            output += "%d\t%d\t%.8e\n" % (src_gid, tgt_gid, w_noise)
                            line_cnt += 1

                            w_mit_pyr_matrix[mit, tgt_mc] = w_noise
                            
        output_fn = self.params['conn_list_mit_pyr']
        print 'Saving %d mit-pyr connections to file %s' % (line_cnt, output_fn)
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(output_fn, 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()

        np.savetxt(self.params['connection_matrix_detailed_ob_oc_dat'], w_mit_pyr_matrix)


    def get_mit_rsnp_connections_from_ALa(self, random_conn=False):
        '''
        This functions creates the OB->OC connections files from the matrix stored in self.params['ob_oc_abstract_weights_fn']
        '''
        print "Drawing MIT - RSNP connections .... "
        abstract_weights = np.loadtxt(self.params['ob_oc_abstract_weights_fn'])
        if random_conn:
            rnd.seed(self.params['random_ob_oc_seed'])
            rnd.shuffle(abstract_weights)
            np.savetxt(self.params['ob_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)
        assert (abstract_weights[:, 0].size == self.params['n_mit'])
        assert (abstract_weights[0, :].size == self.params['n_hc'] * self.params['n_mc'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.min()
        w_mit_rsnp_max = self.params['w_mit_rsnp_max']

        output = ""
        line_cnt = 0
        for tgt_mc in xrange(self.params['n_hc'] * self.params['n_mc']): #column
            for mit in xrange(self.params['n_mit']): # row
                w_in = abstract_weights[mit, tgt_mc]
                if (w_in < 0):
#                    for tgt_rsnp in xrange(int(round(self.params['n_tgt_rsnp_per_mc']))):
#                    tgt_rsnps = utils.get_rnd_targets(self.params['n_rsnp_per_mc'], self.params['n_tgt_rsnp_per_mc'])
                    tgt_rsnps = random.sample(range(self.params['n_rsnp_per_mc']), int(round(self.params['n_tgt_rsnp_per_mc'])))
                    for tgt_rsnp in tgt_rsnps:
                        w_out = (w_in / w_max_abstract) * w_mit_rsnp_max 
                        w_noise = rnd.normal(w_out, w_out * self.params['w_mit_rsnp_sigma_frac'])
#                        w_noise = rnd.normal(w_out, self.params['w_mit_rsnp_max'] * self.params['w_mit_rsnp_sigma_frac'])
                        if (w_noise > self.params['weight_threshold']):
                            src_gid = self.params['mit_offset'] + mit
                            tgt_gid = self.params['rsnp_offset'] + tgt_rsnp + tgt_mc * self.params['n_rsnp_per_mc']
                            output += "%d\t%d\t%.8e\n" % (src_gid, tgt_gid, w_noise)
                            line_cnt += 1

        output_fn = self.params['conn_list_mit_rsnp']
        print 'Saving %d mit-rsnp connections to file %s' % (line_cnt, output_fn)
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(output_fn, 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()

    def get_pyr_readout_connections_from_ALa(self):
        '''
        '''
        print "Drawing OC - Readout connections .... "
        abstract_weights = np.loadtxt(self.params['oc_readout_abstract_weights_fn'])
        assert (abstract_weights[:, 0].size == self.params['n_hc'] * self.params['n_mc'])
        assert (abstract_weights[0, :].size == self.params['n_readout'])
        # scale the abstract weights into the biophysical range
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()
        w_pyr_readout_max = self.params['w_pyr_readout']

        output = ""
        line_cnt = 0
        for tgt_cell in xrange(self.params['n_readout']):
            for src_mc in xrange(self.params['n_hc'] * self.params['n_mc']): # row
                w_in = abstract_weights[src_mc, tgt_cell]
                if (w_in > 0):
                    w_out = (w_in / w_max_abstract) * w_pyr_readout_max
                elif (w_in < 0):
                    w_out = (-1.0) * (w_in / w_min_abstract) * w_pyr_readout_max
                if (abs(w_in > self.params['weight_threshold'])):
                    for src_pyr in xrange(self.params['n_pyr_per_mc']):
                        src_gid = self.params['pyr_offset'] + src_mc * self.params['n_pyr_per_mc'] + src_pyr
                        tgt_gid = self.params['readout_offset'] + tgt_cell
                        output += "%d \t%d\t%.8e\n" % (src_gid, tgt_gid, w_out)
                        line_cnt += 1
        first_line = "%d %d\n" % (line_cnt, 3)
        output_file = open(self.params['conn_list_pyr_readout'], 'w')
        output_file.write(first_line)
        output_file.write(output)
        output_file.close()


    def get_oc_ob_connections_from_ALa(self):
        abstract_weights = np.loadtxt(self.params['oc_ob_abstract_weights_fn'])
        assert (abstract_weights[:,0].shape == self.params['n_hc'] * self.params['n_mc'])
        assert (abstract_weights[0,:].shape == self.params['n_mit'])
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()
        w_pyr_readout_max = self.params['w_pyr_readout']



    def get_oc_oc_connections_from_ALa(self, random_conn=False):
        """
        This is a serial version. To be implemented
        if comm.n_proc > 1: call the parallel version
        else: call the serial version
        """

        print "Drawing OC - OC connections .... "
        abstract_weights = np.loadtxt(self.params['oc_oc_abstract_weights_fn'])
        if random_conn:
            rnd.shuffle(abstract_weights)
            rnd.seed(self.params['random_oc_oc_seed'])
            np.savetxt(self.params['oc_oc_abstract_weights_fn'].rsplit('.dat')[0] + '_random.dat', abstract_weights)

        assert (abstract_weights[:,0].size == self.params['n_hc'] * self.params['n_mc'])
        assert (abstract_weights[0,:].size == self.params['n_hc'] * self.params['n_mc'])
        w_max_abstract = abstract_weights.max()
        w_min_abstract = abstract_weights.min()

        w_pyr_pyr_global_max = self.params['w_pyr_pyr_global_max']
        w_pyr_rsnp_max = self.params['w_pyr_rsnp_max']
        output_pyr_pyr = ""
        line_cnt_pyr_pyr = 0
        output_pyr_rsnp = ""
        line_cnt_pyr_rsnp = 0
        cnt_discarded_conn = 0
        for src_mc in xrange(abstract_weights[:, 0].size):
            for tgt_mc in xrange(abstract_weights[:, 0].size):
                if (src_mc != tgt_mc):
                    w_in = abstract_weights[src_mc, tgt_mc]
                    if (w_in > 0): # draw several pyr -> pyr connections between the two MC
                        src_tgt_dict = {} # src_tgt_dict[src_gid] = [tgt_gid_0, ...] multiple connections between the same source and the same target are forbiddden
                        w_out = (w_in / w_max_abstract) * w_pyr_pyr_global_max
                        src_pyrs = rnd.randint(0, self.params['n_pyr_per_mc'], self.params['n_pyr_pyr_between_2mc'])
                        for src in np.unique(src_pyrs):
                            src_tgt_dict[src] = []
                        for src in src_pyrs:
                            src_pyr = src + src_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                            tgt_pyr = rnd.randint(0, self.params['n_pyr_per_mc']) + tgt_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset']
                            src_tgt_dict[src].append(tgt_pyr)

                        # remove multiple instances of the same src-tgt connection
                        for src in src_pyrs:
                            n1 = len(src_tgt_dict[src])
                            src_tgt_dict[src] = np.unique(src_tgt_dict[src]).tolist()
                            cnt_discarded_conn += n1 - len(src_tgt_dict[src])
                            for tgt_pyr in src_tgt_dict[src]:
                                w_noise = utils.draw_connection(1.0, w_out, noise=self.params['w_pyr_pyr_global_sigma'])
                                if (w_noise > self.params['weight_threshold']):
                                    output_pyr_pyr += "%d %d %.6e\n" % (src_pyr, tgt_pyr, w_noise)
                                    line_cnt_pyr_pyr += 1

                    elif (w_in < 0):
                        w_out = (w_in / w_min_abstract) * w_pyr_rsnp_max
                        src_pyrs = utils.get_rnd_targets(self.params['n_pyr_per_mc'], self.params['n_pyr_rsnp_between_2mc']) # avoid double connections
                        for src in src_pyrs:
                            src_pyr = src + src_mc * self.params['n_pyr_per_mc'] + self.params['pyr_offset'] 
                            tgt_rsnp = rnd.randint(0, self.params['n_rsnp_per_mc']) + tgt_mc * self.params['n_rsnp_per_mc'] + self.params['rsnp_offset']
                            w_noise = utils.draw_connection(1.0, w_out, noise=self.params['w_pyr_rsnp_sigma'])
                            if (w_noise > self.params['weight_threshold']):
                                output_pyr_rsnp += "%d %d %.6e\n" % (src_pyr, tgt_rsnp, w_noise)
                                line_cnt_pyr_rsnp += 1

        print 'Number of discarded pyr-pyr connections:', cnt_discarded_conn
        print 'Number of pyr-rsnp connections:', line_cnt_pyr_rsnp
        print 'Number of pyr-pyr connections:', line_cnt_pyr_pyr
        print 'Number of OC-OC connections:', line_cnt_pyr_pyr + line_cnt_pyr_rsnp
        output_fn_pyr_pyr = self.params['conn_list_pyr_pyr']
        output_file_pyr_pyr = open(output_fn_pyr_pyr, 'w')
        output_file_pyr_pyr.write("%d\t%d\n" % (line_cnt_pyr_pyr, 3))
        output_file_pyr_pyr.write(output_pyr_pyr)
        output_file_pyr_pyr.close()

        output_fn_pyr_rsnp = self.params['conn_list_pyr_rsnp']
        output_file_pyr_rsnp = open(output_fn_pyr_rsnp, 'w')
        output_file_pyr_rsnp.write("%d\t%d\n" % (line_cnt_pyr_rsnp, 3))
        output_file_pyr_rsnp.write(output_pyr_rsnp)
        output_file_pyr_rsnp.close()



    def get_ob_oc_connections_parallel(self):
        '''
        pid 0: Copy the vq_0.csv and the Connection_0_n0.csv files produced by learning algorithm into the data folder.
        Process these two files in order to produce a NEURON readable conn_list output file storing 
        the MIT -> PYR cell connections.
        conn_file: output file from Simon's code, stores the connections strengths in format (row, col) = (tgt_mc, src_mc)
                    to which these minicolumns belong is given by the vq_file
        '''
        if (self.my_rank != 0):
            time.sleep(0.5)
        else:
#            if not (os.path.exists('%s' %  self.params['vq_path'])):
            os.system('cp %s %s' % (self.params['vq_fn'], self.params['vq_path']))
#            if not (os.path.exists('%s' %  self.params['mit_pyr_conn_fn_base'])):
            os.system('cp %s* %s' % (self.params['mit_pyr_conn_fn_base'], self.params['folder_name']))
            self.cat_hebbconn_ob_oc_files() # concatenate the output of the learning algorithms

        # ----- set variable names
        n_mc = self.params['n_mc']
        n_hc = self.params['n_hc']
        n_pyr_per_mc = self.params['n_pyr_per_mc']
        n_rsnp_per_mc = self.params['n_rsnp_per_mc']
        n_orn = self.params['n_orn']
        n_mit = self.params['n_mit']
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        mit_offset = self.params['mit_offset']
        pyr_offset = self.params['pyr_offset']
        rsnp_offset = self.params['rsnp_offset']
        num_connections = 0

        # ------ divide the tgt_mc among the processes: each process deals with a subset of rows in the hebbconnection_ob_oc.txt file
        n_tgt_mcs = n_hc * n_mc         # total number of target MCs
        self.n_my_tgt_mcs = int(n_tgt_mcs / self.n_proc)
        R = n_tgt_mcs % self.n_proc
        offset = min(self.my_rank, R)   # possibly some processes have to deal with more tgt_mcs, i.e. process more lines in the hebbconnection file
        limit_min = int(self.my_rank * self.n_my_tgt_mcs + offset)
        if (self.my_rank < R):
            limit_max = int(limit_min + self.n_my_tgt_mcs + 1)
        else:
            limit_max = int(limit_min + self.n_my_tgt_mcs)
        my_tgt_mcs = range(limit_min, limit_max)

        # ------------- get abstract weights ------------------ 
        abstract_weights_unscaled, abstract_weights = self.get_abstract_weight_matrix_ob_oc()
        w_max = abstract_weights.max()
        w_min = abstract_weights.min()
        self.plot_oc_oc_mat = abstract_weights_unscaled
        # ------------ 'randomize' the connections from OB to OC -------------------------------------------
        """
        According to the output of the abstract VQ and BCPNN learning algorithms, 
        the mitral cells from one glomerulus project all to the same cortical regions.
        This contradicts experimental evidence.
        Thus, we do the following:
        For each glomerulus, the minicolumns in the OC targeted with non-zeros weights in the abstract weight matrix, we choose 

        """
        conn_mat = self.get_vq_result()
        for glom in xrange(self.params['n_mit_y']):
            target_hcs = conn_mat[glom, :].nonzero()
            for tgt_hc in target_hcs:
                for x in xrange(self.params['n_mit_x']):
                    src_index = glom * self.params['n_mit_x'] + x
                    # for each HC choose a random set of minicolumns to which this MT cell projects to by random
                    tgt_mcs = rnd.randint(0, n_mc, n_mc - self.params['n_tgt_mc_per_mit_per_hc']) # can contain e.g. [1, 1, 2] -> effectively only [1, 2]
                    mc_offset = tgt_hc * n_mc
                    for tgt_mc in tgt_mcs:
                        abstract_weights[src_index, tgt_mc + mc_offset] = 0.

        output_fn = self.params['connection_matrix_abstract_ob_oc_dat'].rsplit('.dat')[0] + '_randomized.dat'
        np.savetxt(output_fn, abstract_weights, fmt="%.6E")

        # ------------ draw cell-to-cell connections based on abstract_weights ----------------------------- 
        src_offset = self.params['mit_offset']
        for tgt_mc in my_tgt_mcs:       # process only the ones dedicated to the current process
            tgt_hc = tgt_mc / n_mc
            src_list_output = []        # lists to be written to output file
            tgt_list_output = []
            weight_list_output = []
            n_conn = 0
            for src in xrange(n_mit):
                weight = abstract_weights[src, tgt_mc]
                src_gid = src + src_offset 
                if (weight > 0):
                    # draw a number of target pyr cells within the tgt_mc
                    tgt_list = rnd.randint(0, n_pyr_per_mc, self.params['n_tgt_pyr_per_mc'])
                    tgt_list = np.unique(tgt_list) # avoid multiple connections between 2 cells
                    for tgt in tgt_list:
                        w = (weight / w_max) * self.params['w_mit_pyr_max']
                        w_noise = rnd.normal(w, self.params['w_mit_pyr_max'] * self.params['w_mit_pyr_sigma_frac'])
                        old_sign = np.sign(weight)
                        if (np.sign(w_noise) == old_sign) and  (abs(w_noise) > self.params['weight_threshold']):
                            tgt_offset = tgt_hc * n_mc * n_pyr_per_mc + (tgt_mc % n_mc) * n_pyr_per_mc + pyr_offset
                            tgt_gid = tgt + tgt_offset
                            src_list_output.append(src_gid)
                            tgt_list_output.append(tgt_gid)
                            weight_list_output.append(w_noise)
                            n_conn += 1
                        # else: 
                        #   change in sign -> weight = 0
                        #   weight < weight_treshold -> weight = 0

                elif (weight < 0): # mit -> rsnp
                    # connect to rsnp cells in the tgt_mc
                    tgt_rsnps = rnd.randint(0, self.params['n_rsnp_per_mc'], self.params['n_tgt_rsnp_per_mc'])
                    tgt_rsnps = np.unique(tgt_rsnps) # avoid multiple connections between 2 cells
                    for tgt_rsnp in tgt_rsnps:
                        w = -1. * (weight / w_min) * self.params['w_mit_pyr_min']
                        w_noise = rnd.normal(w, abs(self.params['w_mit_pyr_min']) * self.params['w_mit_pyr_sigma_frac'])
                        old_sign = np.sign(weight)
                        if (np.sign(w_noise) == old_sign) and  (abs(w_noise) > self.params['weight_threshold']):
                            tgt_offset = tgt_hc * n_mc * n_rsnp_per_mc + (tgt_mc % n_mc) * n_rsnp_per_mc + rsnp_offset
                            tgt_gid = tgt_rsnp + tgt_offset
                            src_list_output.append(src_gid)
                            tgt_list_output.append(tgt_gid)
                            weight_list_output.append(abs(w_noise))
                            n_conn += 1
                            if self.debug:
                                assert (src_gid < (self.params['n_orn'] + self.params['n_mit'])), 'Ob->Oc connection wrong, mit gid too high: %d\n' % (src_gid)
                                assert (tgt_gid < (rsnp_offset + self.params['n_rsnp'])), 'Ob->Oc connection wrong, rsnp gid too high: %d\n' % (tgt_gid)
                                assert (tgt_gid >= rsnp_offset), 'Ob->Oc connection wrong, rsnp gid too low: %d\n' % (tgt_gid)

            # write the connections targeting this minicolumn to a file
            output_fn_mc = '%s%d.dat' % (self.params['conn_list_mit_oc_base'], tgt_mc)
            output_data = np.array((src_list_output, tgt_list_output, weight_list_output))
            # binary write
            output_file = open(output_fn_mc,'w')
            output_file.write(output_data.tostring())
            if self.debug == 1: # plain text output
                debug_output_fn = '%s%d_debug.dat' % (self.params['conn_list_mit_oc_base'], tgt_mc)
                np.savetxt(debug_output_fn, output_data.transpose(), fmt='%.10e', delimiter='\t')
            output_file.close()

        if (self.my_rank == 0):
            self.plot_weight_matrix(abstract_weights)
    # ---- end of get_ob_oc_connections_parallel(self) ------------------ 

        
    def plot_weight_matrix(self, matrix):
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print 'plotting ob -> oc weightmatrix as colormap to file: ', self.params['conn_list_mit_oc_colormap']
#        norm = matplotlib.mpl.colors.Normalize(vmin=matrix.min(), vmax=matrix.max())
#        ax.imshow(matrix, interpolation='nearest')
        cax = ax.pcolormesh(matrix)
        pylab.colorbar(cax)
        pylab.savefig(self.params['conn_list_mit_oc_colormap'])
#        pylab.show()
            

    def plot_oc_oc_weight_matrix(self):
        if (self.my_rank == 0):
            fig = pylab.figure()
            ax = fig.add_subplot(111)
            print 'plotting ob -> oc weightmatrix as colormap to file: ', self.params['conn_list_oc_oc_colormap']
            ax.imshow(self.plot_oc_oc_mat, interpolation='nearest')
            pylab.savefig(self.params['conn_list_oc_oc_colormap'])
#            pylab.show()
            


#        self.send_rcv_weight_matrices()
#        if (self.my_rank == 0):
#            fig = pylab.figure()
#            ax = fig.add_subplot(111)
#            print 'plotting ob -> oc weightmatrix as colormap to file: ', self.params['conn_list_mit_oc_colormap']
#            cax = ax.pcolor(self.all_ob_oc_matrices)#, edgecolor='k', linewidths='1')
#            pylab.xlabel('Minicolumn index', fontsize=16)
#            pylab.ylabel('Mitral cell', fontsize=16)
#            pylab.colorbar(cax)
#            pylab.savefig(self.params['conn_list_mit_oc_colormap'])
#            pylab.close()


#        self.all_ob_oc_matrices = None

    def send_rcv_weight_matrices(self):

        n_proc = self.comm.Get_size()
        n_row = self.weight_matrix.shape[0]
        n_col = self.weight_matrix.shape[1]
        all_data = [self.weight_matrix]
        if (self.my_rank == 0):
            for sender in xrange(1, n_proc):
                # each process sends a matrix
                data = self.comm.recv(source=sender, tag=sender)
                all_data.append(data)
                # collect all data into one weight_matrix
                n_row += data.shape[0]
                n_col += data.shape[1]
                
        # all processes except for one sends data
        elif (self.my_rank != 0):
            self.comm.send(self.weight_matrix, dest=0, tag=self.my_rank)

        if (self.my_rank == 0):
            self.all_ob_oc_matrices = np.zeros((n_row, n_col))
            r_old = 0
            c_old = 0
            for p in xrange(len(all_data)): 
                r = all_data[p].shape[0] + r_old # n_rows of this matrix
                c = all_data[p].shape[1] + c_old # n_cols of this matrix
                self.all_ob_oc_matrices[r_old:r, c_old:c] = all_data[p]
                r_old = r
                c_old = c
        all_data = []

    def get_ob_oc_connections_serial(self):
        '''
        Copy the vq_0.csv and the Connection_0_n0.csv files produced by Simon's code into the data folder.
        Process these two files in order to produce a NEURON readable conn_list output file storing 
        the MIT -> PYR cell connections.
        conn_file: output file from Simon's code, stores the connections strengths in format (row, col) = (tgt_mc, src_mc)
                    to which these minicolumns belong is given by the vq_file
        '''
        if (self.my_rank != 0):
            time.sleep(1)
            return
        os.system('cp %s %s' % (self.params['vq_fn'], self.params['vq_path']))
        os.system('cp %s* %s' % (self.params['mit_pyr_conn_fn_base'], self.params['folder_name']))
#        os.system('cp %s %s' % (self.params['mit_pyr_conn_fn'], self.params['mit_pyr_conn_path']))
#        os.system('cp %s %s' % (self.params['pyr_pyr_conn_fn'], self.params['pyr_pyr_conn_path']))

        self.cat_hebbconn_ob_oc_files()

        n_mc = self.params['n_mc']
        n_hc = self.params['n_hc']
        n_pyr_per_mc = self.params['n_pyr_per_mc']
        n_rsnp_per_mc = self.params['n_rsnp_per_mc']
        n_orn = self.params['n_orn']
        n_mit = self.params['n_mit']
        n_mit_x = self.params['n_mit_x']
        n_mit_y = self.params['n_mit_y']
        mit_offset = self.params['mit_offset']
        pyr_offset = self.params['pyr_offset']
        rsnp_offset = self.params['rsnp_offset']

        # set input file names
        conn_file = file(self.params['mit_pyr_conn_path'], 'r')
        print 'Pid %d Creating OB->OC connections from files: %s and %s' % (self.my_rank, self.params['mit_pyr_conn_path'], self.params['vq_path'])
        output_file = open(self.params['conn_list_mit_oc'], 'w')
        output_file_tmp_fn = self.params['conn_list_mit_oc'] + 'tmp'

        vq_file = file(self.params['vq_path'], 'r')
        vq_lines = vq_file.readlines()
        num_tgt_hc = len(vq_lines)
        vq_data = [[] for i in xrange(num_tgt_hc)] 
        weight_matrix = np.zeros((n_mit, n_hc * n_mc))

        magic_number = 123456789
        glom_hc_mapping = np.ones(n_mit_y) # to which tgt_hc does the src_hc connect to?
        glom_hc_mapping *= magic_number
        # hc_hc_mapping[src_hc] = tgt_hc
        for tgt_hc in xrange(num_tgt_hc): #tgt_hc is equivalent to a line in vq_file
            # remove whitespace from the vq lines
            src_hcs = re.sub('\s', '', vq_lines[tgt_hc]) #src_hcs = '2,7'
            if (src_hcs != ''): # if this tgt_hc gets input and is not an empty line in the vq_file ...
                # get the list of source hypercolumns
                src_hcs = src_hcs.rsplit(',') # e.g. src_hcs = [2, 7] 
                for src_hc in src_hcs:
                    glom_hc_mapping[int(src_hc)] = int(tgt_hc)
                    vq_data[tgt_hc].append(int(src_hc))

        # read lines from the Conn_file which stores the weights for
        # tgt_mc vs src_mc (ids for the hcs is stored in vq_file)
        conn_file_lines = conn_file.readlines() # hebbconnections_ob_oc.txt
        mc_mc_conns = [] # a list of lists with conn_file_lines[tgt_mc] = src_mc
        for line_cnt in xrange(len(conn_file_lines)):
            # remove all the whitespace
            conn_file_lines[line_cnt] = re.sub('\s', '', conn_file_lines[line_cnt])
            mc_mc_conns.append(conn_file_lines[line_cnt].rsplit(','))
#        print 'n_lines in conn_file_path:', len(conn_file_lines)

        w_min, w_max = 0.0, 0.0 # just to know what values come out of Simon's code
        for tgt_mc in xrange(len(mc_mc_conns)): # row in hebbconnections 
            for col in xrange(len(mc_mc_conns[tgt_mc])): # col in hebbconnections 
                weight = float(mc_mc_conns[tgt_mc][col])
                if (weight != 0.0):
                    weight = np.log(weight)
                    w_min = min(weight, w_min)
                    w_max = max(weight, w_max)

#        for glom in xrange(n_mit_y):
#            print 'Glom %d -> hc %d' % (mit, glom_hc_mapping[glom])
                # go to the right place in the hebbconnection_ob_oc.txt data

        lines_to_write = '' # output connections
        n_conn = 0
        assert (len(mc_mc_conns) == self.params['n_hc'] * self.params['n_mc']), "Wrong filesize %s\n Clean working directory? no leftovers from earlier (parallel runs?)" % (self.params['mit_pyr_conn_path'])
        for tgt_mc in xrange(len(mc_mc_conns)): # row in hebbconnections 
            tgt_hc = tgt_mc / n_mc
#            if (tgt_mc % n_mc == 0): # create a new folder where all files belonging to this tgt hypercolumn will be stored
#                new_folder = '%s%d' % (self.params['connections_folder_name'], tgt_hc)
#                os.system('mkdir %s' % (new_folder))

            src_list = []
            tgt_list = []
            weight_list = []
            for src_mt in xrange(len(mc_mc_conns[tgt_mc])): # src_mt in hebbconnections 
                # compute the right position in the cell array
#                print 'debug src_mt %d , len(mc_mc_conns[tgt_mc=%d]) = %d, n_mc =%d' % (src_mt, tgt_mc, len(mc_mc_conns[tgt_mc]), n_mc)
#                print 'debug %d tgt_hc %d src_mt / n_mit_x %d , len(vq_data) %d'  % (self.my_rank, tgt_hc, src_mt/n_mit_x, len(vq_data))
#                print 'debug %d len(vq_data[%d]' % (self.my_rank, len(vq_data[tgt_hc]))
#                print 'debug vq_path', self.params['vq_path']
#                print 'debug conn_file', self.params['mit_pyr_conn_path'] 
                glom = vq_data[tgt_hc][src_mt / n_mit_x]
                weight = float(mc_mc_conns[tgt_mc][src_mt])
                if (weight != 0):
                    weight = np.log(weight)
                    # compute the source and target gids
                    src_gid = glom * n_mit_x + src_mt % n_mit_x + mit_offset
                    weight_matrix[src_gid - mit_offset, tgt_hc * n_mc + (tgt_mc % n_mc)] = weight
                    if (weight > 0):
                        # draw a number of target pyr cells within the tgt_mc
                        tgts = rnd.randint(0, n_pyr_per_mc, self.params['n_tgt_pyr_per_mc'])
                        for tgt in tgts:
                            w = (weight / w_max) * self.params['w_mit_pyr_max']
                            w_noise = rnd.normal(w, abs(w) * self.params['w_mit_pyr_sigma_frac'])
                            old_sign = np.sign(weight)
#                            if (np.sign(w_noise) == old_sign):
                            if (np.sign(w_noise) == old_sign) and  (w_noise > self.params['weight_threshold']):
                                tgt_offset = tgt_hc * n_mc * n_pyr_per_mc + (tgt_mc % n_mc) * n_pyr_per_mc + pyr_offset
                                tgt_gid = tgt + tgt_offset
#                                assert (src_gid < (self.params['n_orn'] + self.params['n_mit'])), 'Ob->Oc connection wrong, mit gid too high: %d\n' % (src_gid)
#                                assert (tgt_gid < (pyr_offset + self.params['n_pyr'])), 'Ob->Oc connection wrong, pyr gid too high: %d\n' % (tgt_gid)

                                lines_to_write += '%d\t%d\t%.8e\n' % (src_gid, tgt_gid, w_noise)
                                src_list.append(src_gid)
                                tgt_list.append(tgt_gid)
                                weight_list.append(w_noise)
                                n_conn += 1
                            # else: weight change by adding noise to weight -> no connection drawn
                    elif (weight < 0):
                        # connect to the two rsnp cells in the tgt_mc
                        tgt_rsnp = rnd.randint(0, self.params['n_rsnp_per_mc'])
#                        for tgt_rsnp in xrange(self.params['n_rsnp_per_mc']):
                        w = -1. * (weight / w_min) * self.params['w_mit_pyr_min']
                        w_noise = rnd.normal(w, abs(w) * self.params['w_mit_pyr_sigma_frac'])
                        old_sign = np.sign(weight)
#                        if (np.sign(w_noise) == old_sign):
                        if (np.sign(w_noise) == old_sign) and  (w_noise > self.params['weight_threshold']):
                            tgt_offset = tgt_hc * n_mc * n_rsnp_per_mc + (tgt_mc % n_mc) * n_rsnp_per_mc + rsnp_offset
                            tgt_gid = tgt_rsnp + tgt_offset
                            lines_to_write += '%d\t%d\t%.8e\n' % (src_gid, tgt_gid, w_noise)
                            src_list.append(src_gid)
                            tgt_list.append(tgt_gid)
                            weight_list.append(w_noise)
                            n_conn += 1
#                            assert (src_gid < (self.params['n_orn'] + self.params['n_mit'])), 'Ob->Oc connection wrong, mit gid too high: %d\n' % (src_gid)
#                            assert (tgt_gid < (rsnp_offset + self.params['n_rsnp'])), 'Ob->Oc connection wrong, rsnp gid too high: %d\n' % (tgt_gid)
#                            assert (tgt_gid >= rsnp_offset), 'Ob->Oc connection wrong, rsnp gid too low: %d\n' % (tgt_gid)
                            # else: weight change by adding noise to weight -> no connection drawn



        # first line needed for later processing with NEURON
        first_line = '%d\t3\n' % (n_conn)
        output_file.write(first_line)
        # write the connections targeting this minicolumn to a file
        output_file.write(lines_to_write)
        output_file.close()

        # plot the weight matrix as colormap
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        print 'plotting ob -> oc weightmatrix as colormap to file: ', self.params['conn_list_mit_oc_colormap']
        cax = ax.pcolor(weight_matrix)#, edgecolor='k', linewidths='1')

        pylab.xlabel('Minicolumn index', fontsize=16)
        pylab.ylabel('Mitral cell', fontsize=16)

        pylab.colorbar(cax)
        pylab.savefig(self.params['conn_list_mit_oc_colormap'])
        pylab.close()


     # ---------------- end of get_ob_oc_connections_serial


    def get_oc_readout_connections(self):
        if (self.my_rank != 0):
            time.sleep(1)
            return
        conn_fn = self.params['pyr_readout_conn_path']
        conn_file = file(conn_fn, 'r')
        num_connections = 0
        lines_to_write = '' # output connections

        n_pyr = self.params['n_pyr']
        n_readout = self.params['n_readout']
        pyr_offset = self.params['pyr_offset']
        readout_offset = self.params['readout_offset']
        # assuming n_patterns == n_readout
        # read n_patterns lines from the connection file which stores the weights for
        # the oc -> readout connections
        # line = contains all connections targeting one readout cell
        # line[src_mc] = weight from src_mc to current readout
        conn_file_lines = conn_file.readlines()
        assert (len(conn_file_lines) == self.params['n_readout']), "Wrong filesize: %s\n Cleaned the working directory of other files? " % (conn_fn)
        abstract_weights = np.zeros((n_pyr, n_readout))
        # find out the max, min, mean and median for all pyr - readout weights
        weights = []
        w_max, w_min = 0.0, 0.0
        for tgt_readout in xrange(len(conn_file_lines)):
            conn_file_lines[tgt_readout] = re.sub('\s', '', conn_file_lines[tgt_readout])
            weights = conn_file_lines[tgt_readout].rsplit(',')
            assert (len(weights) == self.params['n_hc'] * self.params['n_mc'])
            for src_mc in xrange(len(weights)):
                w_in = float(weights[src_mc])
                if (w_in != 0.0):
#                if (abs(w_in) > 0.01 * abs(w_max)):
                    w_log = np.log(w_in * self.params['w_pyr_readout_sf'])
                    weights.append(w_log)
                    if (w_log > 0):
                        w_max = max(w_log, w_max)
                    elif (w_log < 0):
                        w_min = min(w_log, w_min)
                    
        for tgt_readout in xrange(len(conn_file_lines)):
            conn_file_lines[tgt_readout] = re.sub('\s', '', conn_file_lines[tgt_readout])
            weights = conn_file_lines[tgt_readout].rsplit(',')
            assert (len(weights) == self.params['n_hc'] * self.params['n_mc'])
            for src_mc in xrange(len(weights)):
                w_in = float(weights[src_mc])
                if (w_in != 0.0):
#                if (abs(w_in) > 0.01 * abs(w_max)):
                    w_log = np.log(w_in * self.params['w_pyr_readout_sf'])
                    if (w_log) > 0:
                        w_out = (w_log / w_max) * self.params['w_pyr_readout']
                    else:
#                        w_out = (-1.0) * (w_log / w_min) * self.params['w_pyr_readout']
                        w_out = 0.0
                    src_offset = src_mc * self.params['n_pyr_per_mc'] + pyr_offset
                    for src in xrange(self.params['n_pyr_per_mc']):
                        src_gid = src + src_offset
                        tgt_gid = tgt_readout + readout_offset
                        lines_to_write += '%d\t%d\t%.8e\n' % (src_gid, tgt_gid, w_out)
                        num_connections += 1
                        
        # first line needed for later processing with NEURON
        first_line = '%d\t3\n' % (num_connections)
        output_file = file(self.params['conn_list_pyr_readout'], 'w')
        output_file.write(first_line)
        output_file.write(lines_to_write)


    def get_abstract_weight_matrix_oc_oc(self):
        """
        Return two connection matrices:
        abstract unscaled 
        abstract scaled and afterwards taking the log
        This function also writes these two connection matrices to the outputfiles:
        """
        fn = self.params['pyr_pyr_conn_path']
        path = os.path.abspath(fn)
        data = np.loadtxt(path, delimiter=",")
        n_row = data[:,0].size
        n_col = data[1,:].size
        weights = np.zeros((n_row, n_col))
        weights_log = np.zeros((n_row, n_col))

        mean_weight = data.mean()
        assert (n_col == self.params['n_hc'] * self.params['n_mc']), "Something went wrong in %s\n Check network params, training_params and if data folder cleaned" % (self.params['pyr_pyr_conn_path'])
        assert (n_row == self.params['n_hc'] * self.params['n_mc']), "Something went wrong in %s\n Check network params, training_params and if data folder cleaned" % (self.params['pyr_pyr_conn_path'])
        weights = np.zeros((self.params['n_hc'] * self.params['n_mc'], self.params['n_hc'] * self.params['n_mc']))
        for i in xrange(n_row):
            for j in xrange(n_col):
                w = data[i, j]
                if (w != 0):
                    weights_log[i, j] = np.log(data[i, j] * self.params['oc_oc_weight_scaling'])
                weights[i, j] = data[i, j]

        output_fn = self.params['connection_matrix_abstract_oc_oc_dat']
        np.savetxt(output_fn, weights_log)
        output_fn_unscaled = self.params['connection_matrix_abstract_oc_oc_dat_unscaled']
        np.savetxt(output_fn_unscaled, weights)
        return weights, weights_log


    def get_oc_oc_connections_parallel(self):
        '''
        This function takes the output of MDS, VQ learning algorithms for the recurrent OC -> OC connections
        and converts the csv files into a NEURON compatible file storing the weights
        for connections between pyramidal neurons in the OC layer.
        '''
        if (self.my_rank != 0):
            time.sleep(1)
        else:
            self.cat_hebbconn_oc_oc_files() # concatenate the output of the learning algorithms

        pyr_offset = self.params['pyr_offset']
        rsnp_offset = self.params['rsnp_offset']
        n_mc = self.params['n_mc']
        n_hc = self.params['n_hc']
        pyr_rsnp_output_fn = self.params['conn_list_pyr_rsnp_base'] + str(self.my_rank) + '.dat'
        pyr_rsnp_output_file = file(pyr_rsnp_output_fn, 'w')

        # ------ divide the tgt_mc among the processes: each process deals with a subset of rows in the hebbconnection_ob_oc.txt file
        n_tgt_mcs = n_hc * n_mc         # total number of target MCs
        self.n_my_tgt_mcs = int(n_tgt_mcs / self.n_proc)
        R = n_tgt_mcs % self.n_proc
        offset = min(self.my_rank, R)   # possibly some processes have to deal with more tgt_mcs, i.e. process more lines in the hebbconnection file
        limit_min = int(self.my_rank * self.n_my_tgt_mcs + offset)
        if (self.my_rank < R):
            limit_max = int(limit_min + self.n_my_tgt_mcs + 1)
        else:
            limit_max = int(limit_min + self.n_my_tgt_mcs)
        my_tgt_mcs = range(limit_min, limit_max)

        # ------- get abstract weights 
        weights, weights_log = self.get_abstract_weight_matrix_oc_oc()
        w_max = weights_log.max()
        w_min = weights_log.min()
        print 'OC -> OC: pyr->pyr w_min : ', w_min, ' w_max: ', w_max
#        assert (w_min < 0), 'Wrong sign for minimal weights between two minicolumns belonging to different patterns: %s' % input_fn

        srcs_pyr_rsnp = []
        tgts_pyr_rsnp = []
        weights_pyr_rsnp = []
        for tgt_mc in xrange(n_mc * n_hc):
            srcs_pyr_pyr = []
            tgts_pyr_pyr = []
            weights_pyr_pyr = []
            row_cnt = 0
            pyr_pyr_output_fn_mc = '%s%d.dat' % (self.params['conn_list_pyr_pyr_base'], tgt_mc)
            pyr_pyr_output_fn_mc_debug = '%s%d_debug.dat' % (self.params['conn_list_pyr_pyr_base'], tgt_mc)
            pyr_pyr_output_file = open(pyr_pyr_output_fn_mc,'wb')
            pyr_pyr_output_file_debug = open(pyr_pyr_output_fn_mc_debug,'w')
            if self.debug:
                pyr_pyr_lines = ""
            for src_mc in xrange(n_mc * n_hc):
                w_in = weights[src_mc, tgt_mc]
                if (src_mc == tgt_mc): # no additional recurrent connections within a minicolumn
                    w_in = 0.0
                if (w_in != 0.0):
                    w_in = np.log(w_in * self.params['oc_oc_weight_scaling']) # = weights_log[src_mc, tgt_mc]
                    if (w_in > 0.0): # pyr->pyr connection
                        # w_out stands for the connection weight between minicolumns
                        # This is now translated into a cell -> cell connection weight
                        w_out = w_in * self.params['w_pyr_pyr_global_max'] / w_max
                        src_offset = pyr_offset + src_mc * self.params['n_pyr_per_mc']
                        tgt_offset = pyr_offset + tgt_mc * self.params['n_pyr_per_mc']
                        n_conn = self.params['n_pyr_pyr_between_2mc']
                        sources = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                        targets = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                        w_with_noise = rnd.normal(w_out, self.params['w_pyr_pyr_global_sigma'], n_conn)

                        for c in xrange(int(n_conn)):
                            src_gid = sources[c] + src_offset
                            tgt_gid = targets[c] + tgt_offset
                            if self.debug:
                                assert ((src_gid >= pyr_offset) and (src_gid < pyr_offset + self.params['n_pyr']))
                                assert ((tgt_gid >= pyr_offset) and (tgt_gid < pyr_offset + self.params['n_pyr']))
                                assert (src_gid < self.params['n_cells'] - self.params['n_readout'])
                                assert (tgt_gid < self.params['n_cells'] - self.params['n_readout'])
                            if ((w_with_noise[c] > self.params['weight_threshold']) and (w_with_noise[c] > 0)):
                                if self.debug:
                                    if (abs(src_gid - tgt_gid) < 30):
                                        print '%d drawing an additional pyr->pyr connection within the same minicolumn %d: %d -> %d, weight=%.2e' % (self.my_rank, tgt_mc, src_gid, tgt_gid, w_with_noise[c])
                                pyr_pyr_lines += '%d\t%d\t%.10e\n' % (src_gid, tgt_gid, w_with_noise[c])
                                srcs_pyr_pyr.append(src_gid)
                                tgts_pyr_pyr.append(tgt_gid)
                                weights_pyr_pyr.append(w_with_noise[c])
                                row_cnt += 1

                    elif (w_in < 0.0): # pyr -> rsnp 
                        assert (src_mc != tgt_mc), 'Trying to draw a pyr -> rsnp connection within the same minicolumn'
                        n_conn = self.params['n_pyr_rsnp_between_2mc'] 
                        sources = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                        targets = rnd.randint(0, self.params['n_rsnp_per_mc'], n_conn)
                        w_out = w_in * self.params['w_pyr_rsnp_max'] / w_min # <0 * >0 / <0
                        w_with_noise = rnd.normal(w_out, self.params['w_pyr_rsnp_sigma'], n_conn)
                        src_offset = pyr_offset + src_mc * self.params['n_pyr_per_mc']
                        tgt_offset = rsnp_offset + tgt_mc * self.params['n_rsnp_per_mc']
                        # w_out should be positive here
                        for c in xrange(int(n_conn)):
#                                if (w_with_noise[c] > 0.0): # it should be positive, otherwise the normal distribution made the change of sign
                            if ((abs(w_with_noise[c]) > self.params['weight_threshold']) and (w_with_noise[c] > 0)):
                                src_gid = sources[c] + src_offset
                                tgt_gid = targets[c] + tgt_offset
                                if self.debug:
                                    assert ((src_gid >= pyr_offset) and (src_gid < pyr_offset + self.params['n_pyr']))
                                    assert (tgt_gid >= pyr_offset + self.params['n_basket'])
                                    assert (src_gid < self.params['n_cells'] - self.params['n_readout'])
                                    assert (tgt_gid < self.params['n_cells'] - self.params['n_readout']), tgt_gid
                                srcs_pyr_rsnp.append(src_gid)
                                tgts_pyr_rsnp.append(tgt_gid)
                                weights_pyr_rsnp.append(w_with_noise[c])
                                # gives 'wrong ' ordering of values in binary file -> remove
                                row_cnt += 1
            # end of for w_in in weights ----------------------------------------

            if self.debug:
                pyr_pyr_output_file_debug.write(pyr_pyr_lines)
                pyr_pyr_output_file_debug.close()
            # after all connection to this MC have been drawn, write the number of lines that have been written to the filesize lookup table
            file_index = self.param_tool.file_counter * self.n_proc + self.my_rank
            self.param_tool.update_file_lookup_tables(pyr_pyr_output_fn_mc, row_cnt, 3, file_index) # write into two files: (filename|index)  and (index|nrow|nrcol)

            output_data = np.array((srcs_pyr_pyr, tgts_pyr_pyr, weights_pyr_pyr))
            pyr_pyr_output_file.write(output_data.tostring())# write the random weights in binary form to the output file
            pyr_pyr_output_file.close()
#            print 'Process %d finished writing into file %s' % (self.my_rank, pyr_pyr_output_fn_mc)
            if (self.my_rank == self.n_proc - 1):
                print 'Get OC->OC connections finished %d percent' % (float(tgt_mc) / (self.params['n_hc'] * self.params['n_mc']) * 100)

#        print 'Process %d writing into file %s' % (self.my_rank, pyr_rsnp_output_fn)
        output_data_pyr_rsnp = np.array((srcs_pyr_rsnp, tgts_pyr_rsnp, weights_pyr_rsnp))
        pyr_rsnp_output_file.write(output_data_pyr_rsnp.tostring())
        pyr_rsnp_output_file.close()

    # ------- end of  get_oc_oc_connections_parallel -------------------- 





    def get_oc_oc_connections_serial(self):
        '''
        This function takes the output of MDS, VQ learning algorithms for the recurrent OC -> OC connections
        and converts the csv files into a NEURON compatible file storing the weights
        for connections between pyramidal neurons in the OC layer.
        '''

        self.cat_hebbconn_oc_oc_files() # concatenate the output of the learning algorithms
        input_fn = self.params['pyr_pyr_conn_path']
        input_file = file(input_fn, 'r')
        print 'Creating OC->OC connections from file: %s' % (input_fn)
        # pyramidal cell gid offset
        pyr_offset = self.params['pyr_offset']
        rsnp_offset = self.params['rsnp_offset']

        # init counters
        pyr_pyr_conns = 0
        pyr_rsnp_conns = 0
        src_mc = 0
        tgt_mc = 0
        pyr_pyr_line = ''
        pyr_rsnp_line = ''

        w_max = 0.
        w_min = 0.
        # iterate through all lines to find w_max and w_min
        for line in input_file:
            weights = line.rsplit(',')
            for w_in in weights: # w_in is the one given by Simon's learning output
#                print 'pid ', self.my_rank, 'w_in:', w_in
                if (float(w_in) != 0):
                    w_in = np.log(float(w_in))
                    w_max = max(w_in, w_max)
                    w_min = min(w_in, w_min)
        input_file.close()
        input_file = file(input_fn, 'r')
        print 'OC -> OC: pyr->pyr w_min : ', w_min, ' w_max: ', w_max
#        assert (w_min < 0), 'Wrong sign for minimal weights between two minicolumns belonging to different patterns: %s' % input_fn

        pyr_rsnp_output_fn = self.params['conn_list_pyr_rsnp_base'] + str(self.my_rank) + '.dat'
        pyr_rsnp_output_file = file(pyr_rsnp_output_fn, 'w')
        pyr_pyr_output_fn = self.params['conn_list_pyr_pyr']
        pyr_pyr_output_file = open(pyr_pyr_output_fn,'w')
        lines_to_write_pyr_pyr = ''
        lines_to_write_pyr_rsnp = ''
        # iterate through all lines
        for line in input_file: # is equivalent to a tgt_mc in cortex
            src_mc = 0
            if ((tgt_mc % self.params['n_proc']) == self.my_rank):
                row_cnt = 0
#                if (self.my_rank == 1):
                print 'Pid: %d, line_cnt or tgt_mc: %d, Percent finish: %d'   % (self.my_rank, tgt_mc, 100. * tgt_mc / (self.params['n_mc'] * self.params['n_hc']))

                # output_filename for pyr - pyr connections
                weights = line.rsplit(',') # convert string of weights in that line into a list
                for w_in in weights: # w_in is the one given by Simon's learning output
                    w_in = float(w_in)
                    if (float(w_in) != 0.0):
    #                    print 'debug, src_mc,', src_mc, 'tgt_mc', tgt_mc, 'w_in', w_in, 'log(win)', np.log(w_in)
                        w_in = np.log(w_in)
                        if (w_in > 0.0): # pyr->pyr connection
                            # w_out stands for the connection weight between minicolumns
                            # This is now translated into a cell -> cell connection weight
                            w_out = w_in * self.params['w_pyr_pyr_global_max'] / w_max
                            src_offset = pyr_offset + src_mc * self.params['n_pyr_per_mc']
                            tgt_offset = pyr_offset + tgt_mc * self.params['n_pyr_per_mc']
                            n_conn = self.params['n_pyr_pyr_between_2mc']
                            sources = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                            targets = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                            w_with_noise = rnd.normal(w_out, self.params['w_pyr_pyr_global_sigma'], n_conn)

                            for c in xrange(int(n_conn)):
                                src_gid = sources[c] + src_offset
                                tgt_gid = targets[c] + tgt_offset
#                                assert ((src_gid >= pyr_offset) and (src_gid < pyr_offset + self.params['n_pyr']))
#                                assert ((tgt_gid >= pyr_offset) and (tgt_gid < pyr_offset + self.params['n_pyr']))
#                                assert (src_gid < self.params['n_cells'] - self.params['n_readout'])
#                                assert (tgt_gid < self.params['n_cells'] - self.params['n_readout'])
#                                if (w_with_noise[c] > 0.0): 
                                if (w_with_noise[c] > self.params['weight_threshold']):
#                                    if (abs(src_gid - tgt_gid) < 30):
#                                        print '%d drawing an additional pyr->pyr connection within the same minicolumn: %d -> %d, weight=%.2e' % (self.my_rank, src_gid, tgt_gid, w_with_noise)
                                    lines_to_write_pyr_pyr += '%d\t%d\t%.10e\n' % (src_gid, tgt_gid, w_with_noise[c])
                                    pyr_pyr_conns += 1
                        if (w_in < 0.0): # pyr -> rsnp 
                            assert (src_mc != tgt_mc), 'Trying to draw a pyr -> rsnp connection within the same minicolumn'
                            n_conn = self.params['n_pyr_rsnp_between_2mc'] 
                            sources = rnd.randint(0, self.params['n_pyr_per_mc'], n_conn)
                            targets = rnd.randint(0, self.params['n_rsnp_per_mc'], n_conn)
                            w_out = w_in * self.params['w_pyr_rsnp_max'] / w_min # <0 * >0 / <0
                            w_with_noise = rnd.normal(w_out, self.params['w_pyr_rsnp_sigma'], n_conn)
                            src_offset = pyr_offset + src_mc * self.params['n_pyr_per_mc']
                            tgt_offset = rsnp_offset + tgt_mc * self.params['n_rsnp_per_mc']
                            # w_out should be positive here
                            for c in xrange(int(n_conn)):
#                                if (w_with_noise[c] > 0.0): # it should be positive, otherwise the normal distribution made the change of sign
                                if (w_with_noise[c] > self.params['weight_threshold']):
                                    src_gid = sources[c] + src_offset
                                    tgt_gid = targets[c] + tgt_offset
#                                    assert ((src_gid >= pyr_offset) and (src_gid < pyr_offset + self.params['n_pyr']))
#                                    assert (tgt_gid >= pyr_offset + self.params['n_basket'])
#                                    assert (src_gid < self.params['n_cells'] - self.params['n_readout'])
#                                    assert (tgt_gid < self.params['n_cells'] - self.params['n_readout']), tgt_gid
#                                    pyr_rsnp_line = '%d\t%d\t%.10e\n' % (src_gid, tgt_gid, w_with_noise[c])
#                                    pyr_rsnp_output_file_tmp.write(pyr_rsnp_line)
                                    lines_to_write_pyr_rsnp += '%d\t%d\t%.10e\n' % (src_gid, tgt_gid, w_with_noise[c])
                                    pyr_rsnp_conns += 1
    #                else:
    #                    print 'zero weight from mc ', src_mc, ' to mc', tgt_mc
                    src_mc += 1
                tgt_mc += 1
                # after all connection to this MC have been drawn, write the number of lines that have been written to the filesize lookup table

                file_index = self.param_tool.file_counter * (self.my_rank + 1)
                self.param_tool.update_file_lookup_tables(pyr_pyr_output_fn_mc, row_cnt, 3, file_index) # write into two files: (filename|index)  and (index|nrow|nrcol)
            else:
                tgt_mc += 1
#            tgt_mc = tgt_mc % self.params['n_mc'] # this is wrong, as mcs are counted continuously and NOT between 0 .. n_mc , thus tgt_mc can be larger than n_mc


        print 'Finished writing into file %s' % (pyr_rsnp_output_fn)
        print 'Finished writing into file %s' % (pyr_pyr_output_fn)
        pyr_pyr_output_file.write('%d\t%d\n' % (pyr_pyr_conns, 3))
        pyr_rsnp_output_file.write('%d\t%d\n' % (pyr_rsnp_conns, 3))
        pyr_pyr_output_file.write(lines_to_write_pyr_pyr)
        pyr_rsnp_output_file.write(lines_to_write_pyr_rsnp)
        pyr_pyr_output_file.close()
        pyr_rsnp_output_file.close()
        


#        pyr_rsnp_output_file_tmp.close()
#        pyr_rsnp_first_line = '%d %d\n' % (pyr_rsnp_conns, 3)
#        pyr_rsnp_output_file.write(pyr_rsnp_first_line)
#        pyr_rsnp_output_file_tmp = file(pyr_rsnp_output_fn_tmp, 'r')
#        written_lines = pyr_rsnp_output_file_tmp.readlines()
#        for line in written_lines:
#            pyr_rsnp_output_file.write(line)
#        pyr_rsnp_output_file.close()
#        pyr_rsnp_output_file_tmp.close()
#        os.system('rm %s' % (pyr_rsnp_output_fn_tmp))


    # ------- end of  get_oc_oc_connections_serial -------------------- 



    def merge_conn_file(self, basename, output_fn):
        '''
        basename : filename basis equal for all processes
        merge the connection files which have been written by several processes into one
        for each file: (written by an individual process)
        1) read first line, increase total number of connections to written in the final output file
        2) write_all_other files to a temp
        '''
        n_conn_total = 0
        if (self.my_rank != 0):
            time.sleep(1)
            return

        # read first line from each of the tmp files and keep only the body in another tmp file
        for p in xrange(self.params['n_proc']):
            input_filename = basename + str(p) + '.dat'
            tmp_fn = input_filename + 'tmp'
            output_file_tmp = open(tmp_fn, 'w')
            print 'Merging conn files ... ', input_filename 
            f = open(input_filename, 'r')

            all_lines = f.readlines()
            # remove whitespace
            first_line = re.sub('\s', '', all_lines[0]) # first_line = '1000, 3'
            n_conn = first_line[0]
            n_col = first_line[1]
            n_conn_total += int(n_conn)
            for l in xrange(1, len(all_lines)):
                output_file_tmp.write(all_lines[l])

            output_file_tmp.close()
        
        # write the first line into final output file
        output_file = open(output_fn, 'w')
        output_file.write('%d %d\n' % (n_conn_total, int(n_col)))
        output_file.close()
        cat_command = 'cat '

        for p in xrange(self.params['n_proc']):
            input_filename = basename + str(p) + '.dat'
            tmp_fn = input_filename + 'tmp'
            cat_command += tmp_fn + ' '

        cat_command += ' > %s'  % (output_fn)
        print 'Pid %d: %s' % (self.my_rank, cat_command)
        os.system(cat_command)

        for p in xrange(self.params['n_proc']):
            input_filename = basename + str(p) + '.dat'
            tmp_fn = input_filename + 'tmp'
            os.system('rm %s' % tmp_fn) 


    def get_matrix_from_conn_list(self, conn_list_fn, src_type, tgt_type):

        if ((src_type == 'pyr') or (src_type == 'rsnp')):
            n_src = self.params['n_mc'] * self.params['n_hc']
        else:
            n_src = self.params['n_%s' % src_type]
        if ((tgt_type == 'pyr') or (tgt_type == 'rsnp')):
            n_tgt = self.params['n_mc'] * self.params['n_hc']
        else:
            n_tgt = self.params['n_%s' % tgt_type]
        src_offset = self.params['%s_offset' % src_type]
        tgt_offset = self.params['%s_offset' % tgt_type]

        m = np.zeros((n_src, n_tgt))
        d = np.loadtxt(conn_list_fn, skiprows=1)

        if ((tgt_type == 'pyr') or (tgt_type == 'rsnp')):
            n_tgt_per_mc = self.params['n_%s_per_mc' % tgt_type]
        else:
            n_tgt_per_mc = 1
        if ((src_type == 'pyr') or (src_type == 'rsnp')):
            n_src_per_mc = self.params['n_%s_per_mc' % src_type]
        else:
            n_src_per_mc = 1

        for i in xrange(d[:, 0].size):
            src = (d[i, 0] - src_offset) / n_src_per_mc
            tgt = (d[i, 1] - tgt_offset) / n_tgt_per_mc
            m[src, tgt] += d[i, 2]

        return m
