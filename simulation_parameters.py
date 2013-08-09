import json
import numpy as np
import numpy.random as rnd
import os

class parameter_storage(object):
    """
    This class contains the simulation parameters in a dictionary called params.
    """

    def __init__(self):

        self.params = {}
        self.set_default_params()
        self.set_filenames()

    def set_default_params(self):


        self.params['OR_activation_normalization'] = False
        self.params['with_artificial_orns'] = 0
        
        self.params['Cluster'] = 0
        self.params['concentration_sweep'] = 0
        self.params['n_patterns'] = 10
        self.params['n_proc'] = 8   # on how many processors do you want to run the neuron code?

        self.params['with_noise'] = 1

        # ------ S E E D S -------------
        self.params['seed_activation_matrix'] = 312
        self.params['seed'] = 98765 # this is for pattern generation, weight randomization, etc
        self.params['netstim_seed'] = 0 # netstim_seed acts as an offset for the RNG provided to the noise input via NetStims

        # ------ N E T W O R K    S I Z E -----------
        if (self.params['concentration_sweep'] == 1):  # it's a concentration sweep
            self.params['n_patterns'] = 1
            if (self.params['Cluster'] == 1):
                self.params['n_or'] = 64
            else:
                self.params['n_or'] = 16
        else:
#            self.params['n_or'] = 16
            self.params['n_or'] = 40
#            self.params['n_or'] = self.params['n_patterns']
        if (self.params['Cluster'] == 1):
            self.params['rel_orn_mit'] = 200
            self.params['rel_gran_mit'] = 100# number of granule cells per mitral cell
            self.params['rel_pg_mit']  = 20# number of periglomerular cells per mitral cell, ~ 20 according to Shepherd
        else:
            self.params['rel_orn_mit'] = 5
            self.params['rel_gran_mit'] = 10# number of granule cells per mitral cell
            self.params['rel_pg_mit']  = 10# number of periglomerular cells per mitral cell, ~ 20 according to Shepherd

        # ------ C E L L     N U M B E R S ------------
        self.params['n_gor'] = 10# number of mitral cells per glomerulus
#        self.params['n_gor'] = 8# number of mitral cells per glomerulus
        self.params['n_glom'] = self.params['n_or']
        self.params['n_orn_x'] = self.params['n_gor'] * self.params['rel_orn_mit']# n_orn_x : number of orns expressing one olfactory receptor
        self.params['n_orn_y'] = self.params['n_or']# n_orn_y : number of different receptor families (each having a different affinity to an odour)
        self.params['n_mit_x'] = self.params['n_gor']	# n_mit_x : number of mitral cells per glomerulus
        self.params['n_mit_y'] = self.params['n_or']	# n_mit_y : number of glomeruli, or hypercolumns on the OB level
        self.params['n_pg_x']  = self.params['rel_pg_mit'] * self.params['n_mit_x']
        # PG cells are divided in several different subpopulations, inspiration came from Toida2008: first distinction is in PG cells making only serial or only reciprocal synapses
        self.params['rel_serial_reciprocal'] = 0.5 # each mitral cell has a bunch of pg cells surrounding it, 'rel_serial_reciprocal' determines how many pg cells make serial and how many reciprocal synapses
        self.params['rel_reciprocal_intraglom'] = 0.1 # most of the PG cells in the reciprocal group make only local connections, but 'rel_reciprocal_intraglom' of this group do make ddi connections also with other mitral cells in the glomerulus
        self.params['rel_ff_wta_inhibition'] = 0.7 # rel_ff_wta_inhibition * n_pg_x_serial get excitation from one mitral cell but inhibit all other mitral cells via serial connections
        # (1 - rel_ff_wta_inhibition) * n_pg_x_serial get excitation from ORNs and provide feed-forward inhibition
        self.params['n_pg_x_serial'] = int(round(self.params['n_pg_x'] * self.params['rel_serial_reciprocal']))
        self.params['n_pg_x_rec'] = self.params['n_pg_x'] - self.params['n_pg_x_serial']
        #constrained by n_mit_x and prop_pg_mit_serial_rec, detemined numerically from a plot to get the desired prop_pg_mit_serial_rec
        # checkout ~/workspace/log/111103_OB_connectivity_revisited, Shepherd Synaptic Organization of the brain says: ca. 25% of the dendro-dendritic synapses in the glomerular layer are reciprocal
        self.params['n_pg_y']  = self.params['n_glom']
        self.params['n_gran_x'] = self.params['rel_gran_mit'] * self.params['n_mit_x'] 
        self.params['n_gran_y'] = self.params['n_glom']
        self.params['n_orn'] = self.params['n_orn_x'] * self.params['n_orn_y']
        self.params['n_mit'] = self.params['n_mit_x'] * self.params['n_mit_y']
        self.params['n_pg']  = self.params['n_pg_x'] * self.params['n_pg_y']
        self.params['n_gran'] = self.params['n_gran_x'] * self.params['n_gran_y']
        self.params['n_cells_ob'] = self.params['n_mit'] + self.params['n_gran'] + self.params['n_pg']
        self.params['prop_pg_mit_serial_rec'] = 3 # relation between number serial and reciprocal synapses between periglomerular and MT cells, 3 is according to Shepherd's Book "Synaptic Organization of the Brain",i.e. 25% are reciprocal synapses

        self.params['n_hc'] = 20
        self.params['n_mc'] = 20
        self.params['n_tgt_basket_per_mc'] = 8 # pyr within one minicolumn connect to this number of 'closest' basket cells
        self.params['n_basket_per_mc'] = 6 #this does not mean that the basket cell is exclusively for the minicolumn
        self.params['n_basket_per_hc'] = self.params['n_mc'] * self.params['n_basket_per_mc']
        self.params['n_pyr_per_mc'] = 30
#        self.params['n_tgt_mc_per_mit_per_hc'] = int(round(self.params['n_mc'] / 4.))
        self.params['n_tgt_pyr_per_mc'] = self.params['n_pyr_per_mc'] / 2.0 # number of pyr cells per minicolumn activated by input from OB
#        self.params['n_pyr_pyr_between_2mc'] =  self.params['n_hc'] * self.params['n_pyr_per_mc'] * 0.33 # number of pyr->pyr connections between two minicolumns (belonging to the same pattern
        self.params['n_pyr_pyr_between_2mc'] =  self.params['n_pyr_per_mc'] ** 2 * 0.05 # number of pyr->pyr connections between two minicolumns (belonging to the same pattern)
        self.params['n_pyr_rsnp_between_2mc'] = self.params['n_pyr_per_mc'] / 3.0 # number of pyr->rsnp connections between two minicolumns (belonging to different patterns)
        self.params['n_rsnp_per_mc'] = 4
        self.params['n_tgt_rsnp_per_mc'] = self.params['n_rsnp_per_mc'] * .75 # number of MT -> rsnp cell connections if MT cell has inhibitory connection to that minicolumn
        self.params['n_pyr'] = self.params['n_mc'] * self.params['n_hc'] * self.params['n_pyr_per_mc']
        self.params['n_basket'] = self.params['n_hc'] * self.params['n_basket_per_hc']
        self.params['n_rsnp'] = self.params['n_mc'] * self.params['n_hc'] * self.params['n_rsnp_per_mc']
        self.params['n_readout'] = self.params['n_patterns']
        self.params['n_cells_oc'] = self.params['n_pyr'] + self.params['n_rsnp'] + self.params['n_basket'] + self.params['n_readout']
        # gid offsets for various cell types
        self.params['global_offset'] = 0 # GID start value
        self.params['orn_offset'] = self.params['global_offset']
        self.params['mit_offset'] = self.params['orn_offset'] + self.params['n_orn']
        self.params['gran_offset'] = self.params['mit_offset'] + self.params['n_mit']
        self.params['pg_offset'] = self.params['gran_offset'] + self.params['n_gran']
        self.params['pyr_offset'] = self.params['pg_offset'] + self.params['n_pg']
        self.params['basket_offset'] = self.params['pyr_offset'] + self.params['n_pyr']
        self.params['rsnp_offset'] = self.params['basket_offset'] + self.params['n_basket']
        self.params['readout_offset'] = self.params['rsnp_offset'] + self.params['n_rsnp']
        self.params['cell_types'] = ['orn', 'mit', 'gran', 'pg', 'pyr', 'basket', 'rsnp', 'readout']

        # cell gid ranges for each celltype 
        self.params['orn_range'] = (self.params['orn_offset'], self.params['orn_offset'] + self.params['n_orn'])
        self.params['mit_range'] = (self.params['mit_offset'], self.params['mit_offset'] + self.params['n_mit'])
        self.params['gran_range'] = (self.params['gran_offset'], self.params['gran_offset'] + self.params['n_gran'])
        self.params['pg_range'] = (self.params['pg_offset'], self.params['pg_offset'] + self.params['n_pg'])
        self.params['pyr_range'] = (self.params['pyr_offset'], self.params['pyr_offset'] + self.params['n_pyr'])
        self.params['basket_range'] = (self.params['basket_offset'], self.params['basket_offset'] + self.params['n_basket'])
        self.params['rsnp_range'] = (self.params['rsnp_offset'], self.params['rsnp_offset'] + self.params['n_rsnp'])
        self.params['readout_range'] = (self.params['readout_offset'], self.params['readout_offset'] + self.params['n_readout'])


        # n_cells : total number of cells
        self.params['n_cells'] =  self.params['n_orn'] + self.params['n_mit'] + self.params['n_pg'] + self.params['n_gran'] + \
                self.params['n_pyr'] + self.params['n_basket'] + self.params['n_rsnp'] + self.params['n_readout']
        self.params['n_cell_per_glom'] = self.params['n_orn_x'] + self.params['n_mit_x'] + self.params['n_pg_x'] + self.params['n_gran_x']

        # number of randomly selected testcells from which membrane potentials will be recorded
        self.params['n_test_orn'] = 5
        self.params['n_test_mit'] = 5
        self.params['n_test_gran'] = 5
        self.params['n_test_pg'] = 5
        self.params['n_test_pyr'] = 5
        self.params['n_sample_pyr_per_mc'] = 1
        self.params['n_test_basket'] = 1
        self.params['n_sample_basket_per_hc'] = 1
#        self.params['n_sample_basket_per_hc'] = int(round(self.params['n_basket_per_hc'] / 10.0))
        self.params['n_test_rsnp'] = 1
        self.params['n_sample_rsnp_per_mc'] = 1



        # ODORANT - OR DISTRIBUTION PARAMERS
        # obtained through average_OR_affinity_distributions.py
        # The odorant_receptor_distance_range marks the range of possible distances between ORs and odorants based
        # on the clustering results obtained from average_OR_affinity_distributions.py
        self.params['odorant_receptor_distance_range'] = (0, 4.330310991999920844e+01)
        self.params['odorant_receptor_distance_distribution_parameters'] = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
         1.909487e+01, 1.110809e+01, 3.353855e+00, \
         4.188897e+00, 4.088460e+01, 4.966478e-01] # these values are taken from clustering the odorant space with 40 ORs

        """
        # ---------------- C E L L    P A R A M E T E R S --------- # 
        # ---------------- ORN cell parameters:
        """

        # gor stands for the maximum conductance evoked by an odor
        # gor values are distributed between a min and max value
        # good values for system without noise
        self.params['gor_min'] = 3e-5 
        self.params['gor_max'] = 5e-4
        self.params['gor_exp'] = 3
        # if ORNs have all the same conductance parameters, this is the list:
        self.params['gna'] = 0.5        # [S/cm2]
        self.params['gk'] = 0.05        # [S/cm2]
        self.params['gkcag'] = 0.01     # [S/cm2]
        self.params['gcal'] = 6e-4      # [S/cm2]
        self.params['gleak_orn'] = 1e-4 # [S/cm2]
        self.params['tau_cadec'] = 1000 # [ms]

        # parameters for gleak, gkcag, gcal, gained through combined hand tuning / fitting procedure 
#        self.params['gkcag_params'] = [4.99086530e-03, 2.26738160e-02, 2.26738160e-02]
#        self.params['gcal_params'] =  [4.99086531e-04, 2.26738160e-03, 2.26738160e-03]
#        self.params['gleak_params'] = [4.25453912e-05, -5.18713818e+05, 3.47077557e-05]

        self.params['gkcag_params'] = [4.99086530e-03, 2.26738160e-02, 2.26738160e-02]
#        self.params['gcal_params'] =  [4.99086531e-04, 2.26738160e-03, 2.26738160e-03]
        self.params['gcal_params'] =  [4.99086531e-04, 2.0e-3]
        self.params['gleak_params'] = [5.0e-05, -6.0e+05, 1.0e-05]


        # ---------------- E X P E R I M E N T A L    P A R A M E T E R S --------- #
        self.params['temperature'] = 36# [Celsius] very important! required for NEURON simulations
        self.params['t_sim']	= 1200 # [ms] simulated time
        self.params['time_step']= 0.025   # [ms] max time step
        self.params['time_step_rec']= 0.5  # [ms] time step for recording membrane potentials etc
        self.params['thresh']	= 0     # [mV] threshold for spike detection. thresh is currently the same for all cells 	
        self.params['t_start']	= 0     # [ms] start time of current injection into orn cells
        self.params['tau_odorinput_sigmoid'] = 20 # [ms] time constant for sigmoidal function for input conductance time course, check with neuron_files/odorinput.mod
        self.params['t_stop']	= 25 * self.params['tau_odorinput_sigmoid'] # [ms] start time for decaying sigmoid for odor input conductance
#        self.params['curr_amp']= 100	# [nA] amplitude of current injected into orn cells
        self.params['v_init'] = -70.	# [mV]
        self.params['v_init_sigma'] = 3 # [mV]

        # ---------------- N O I S E    P A R A M E T E R S ----------------- # 
        # dummy noise parameters
#        self.params['f_exc_noise_orn'] = 1e-6  # [Hz] Noise inserted into ORN soma generated by NEURON's NetStim mechanism, i.e. Poisson noise in this case
#        self.params['f_inh_noise_orn'] = 1e-6  # [Hz] Poisson input spike trains can be inserted into cells
#        self.params['w_exc_noise_orn'] = 1e-12
#        self.params['w_inh_noise_orn'] = 1e-12
#        self.params['f_exc_noise_mit'] = 1e-6
#        self.params['f_inh_noise_mit'] = 1e-6
#        self.params['w_exc_noise_mit'] = 1e-12
#        self.params['w_inh_noise_mit'] = 1e-12
#        self.params['f_exc_noise_gran'] = 1e-6
#        self.params['f_inh_noise_gran'] = 1e-6
#        self.params['w_exc_noise_gran'] = 1e-12
#        self.params['w_inh_noise_gran'] = 1e-12
#        self.params['f_exc_noise_pg'] = 1e-6
#        self.params['f_inh_noise_pg'] = 1e-6
#        self.params['w_exc_noise_pg'] = 1e-12
#        self.params['w_inh_noise_pg'] = 1e-12
#        self.params['f_exc_noise_pyr'] = 1e-6
#        self.params['f_inh_noise_pyr'] = 1e-6
#        self.params['w_exc_noise_pyr'] = 1e-12
#        self.params['w_inh_noise_pyr'] = 1e-12
#        self.params['f_exc_noise_basket'] = 1e-6
#        self.params['f_inh_noise_basket'] = 1e-6
#        self.params['w_exc_noise_basket'] = 1e-12
#        self.params['w_inh_noise_basket'] = 1e-12
#        self.params['f_exc_noise_rsnp'] = 1e-6
#        self.params['f_inh_noise_rsnp'] = 1e-6
#        self.params['w_exc_noise_rsnp'] = 1e-12
#        self.params['w_inh_noise_rsnp'] = 1e-12


        # real noise parameters
        self.params['f_exc_noise_orn'] = 400.  # [Hz] Noise inserted into ORN soma generated by NEURON's NetStim mechanism, i.e. Poisson noise in this case
        self.params['f_inh_noise_orn'] = 400.  # [Hz] Poisson input spike trains can be inserted into cells
        self.params['w_exc_noise_orn'] = 0.0005 # [uS] exc noise targetting ORN soma
        self.params['w_inh_noise_orn'] = 0.003  # [uS] inh noise targetting ORN soma
        self.params['f_exc_noise_mit'] = 400.
        self.params['f_inh_noise_mit'] = 400.
        self.params['w_exc_noise_mit'] = 0.0005
        self.params['w_inh_noise_mit'] = 0.003 
        self.params['f_exc_noise_gran'] = 400.
        self.params['f_inh_noise_gran'] = 400.
        self.params['w_exc_noise_gran'] = 0.0005
        self.params['w_inh_noise_gran'] = 0.003
        self.params['f_exc_noise_pg'] = 400.
        self.params['f_inh_noise_pg'] = 400.
        self.params['w_exc_noise_pg'] = 0.0005
        self.params['w_inh_noise_pg'] = 0.003
        self.params['f_exc_noise_pyr'] = 400.
        self.params['f_inh_noise_pyr'] = 400.
        self.params['w_exc_noise_pyr'] = 0.002
        self.params['w_inh_noise_pyr'] = 0.001
        self.params['f_exc_noise_basket'] = 400.
        self.params['f_inh_noise_basket'] = 400.
        self.params['w_exc_noise_basket'] = 0.002
        self.params['w_inh_noise_basket'] = 0.001
        self.params['f_exc_noise_rsnp'] = 400.
        self.params['f_inh_noise_rsnp'] = 400.
        self.params['w_exc_noise_rsnp'] = 0.002
        self.params['w_inh_noise_rsnp'] = 0.001



    def set_folder_name(self, folder_name=None):

        folder_name = 'HandTuned'
        self.params['folder_name'] = os.path.abspath(folder_name)
        print 'Folder name:', self.params['folder_name']


    def set_filenames(self, folder_name=None):

        self.set_folder_name(folder_name)
        print 'Folder name:', self.params['folder_name']

        # FOLDER NAMES
        self.params['conn_folder'] = '%s/Connections' % (self.params['folder_name']) # for conn_lists_ only
        self.params['params_folder'] = '%s/Parameters' % (self.params['folder_name']) # for cell parameters only
        self.params['spiketimes_folder'] = '%s/Spiketimes' % (self.params['folder_name']) # for spiketimes
        self.params['nspikes_folder'] = '%s/NumberOfSpikes' % (self.params['folder_name']) # for files storing only the number of spikes
        self.params['volt_folder'] = '%s/VoltageTraces' % (self.params['folder_name']) # for voltage and calcium traces
        self.params['isyn_folder'] = '%s/SynapticCurrents' % (self.params['folder_name']) # for voltage and calcium traces
        self.params['netcon_folder'] = '%s/Netcons' % (self.params['folder_name']) # for 'offline simulations' storing input files with netcon, weight, spiketime etc
        self.params['bcpnn_folder'] = '%s/Bcpnn' % (self.params['folder_name']) # for output of bcpnn algorithm
        self.params['figure_folder'] = '%s/Figures' % (self.params['folder_name']) # for other stuff
        self.params['other_folder'] = '%s/Other' % (self.params['folder_name']) # for other stuff
        self.params['tmp_folder'] = '%s/Tmp' % (self.params['folder_name']) # for tmp stuff
        self.params['input_spikes_folder'] = '%s/InputSpiketrains' % (self.params['folder_name']) # for tmp stuff
        self.params['folder_names'] = [self.params['conn_folder'], \
                                        self.params['params_folder'], \
                                        self.params['spiketimes_folder'], \
                                        self.params['nspikes_folder'], \
                                        self.params['volt_folder'], \
                                        self.params['isyn_folder'], \
                                        self.params['netcon_folder'], \
                                        self.params['bcpnn_folder'], \
                                        self.params['figure_folder'], \
                                        self.params['tmp_folder'], \
                                        self.params['input_spikes_folder'], \
                                        self.params['other_folder']]
        self.create_folders()

        self.params['params_fn_json'] = '%s/simulation_parameters.json' % (self.params['params_folder'])
        self.params['hoc_file'] = '%s/simulation_params.hoc' % (os.path.abspath(self.params['params_folder']))
        self.params['info_file'] =  '%s/simulation_parameters.info' % (self.params['folder_name']) # human readable file, format: parameter = value

        self.params['activation_matrix_fn'] = '%s/activation_matrix.dat' % (self.params['params_folder'])
        self.params['activation_matrix_fig'] = '%s/activation_matrix.png' % (self.params['figure_folder'])

        # parameter files
        self.params['orn_params_fn_base'] =  '%s/orn_params_' % (self.params['params_folder'])
        self.params['mit_params_fn_base'] =  '%s/mit_params_' % ( self.params['params_folder'])
        self.params['pyr_params_file'] =  '%s/pyr_params.dat' % ( self.params['params_folder'])
        self.params['basket_params_fn_base'] =  '%s/basket_params_' % ( self.params['params_folder'])
        self.params['rsnp_params_fn_base'] =  '%s/rsnp_params_' % ( self.params['params_folder'])
        self.params['readout_params_file'] =  '%s/readout_params.dat' % ( self.params['params_folder'])


        # connectivity files
        self.params['conn_list_orn_mit'] =  '%s/conn_list_orn_mit.dat' % ( self.params['conn_folder'])
        self.params['conn_list_orn_pg'] =  '%s/conn_list_orn_pg.dat' % ( self.params['conn_folder']) # for normal connectivity


        # files for recording currents, membrane potential, time, ....
        self.params['time_vec_fn_base'] = '%s/time_vector' % ( self.params['volt_folder'])
        # voltage files: each cell has its own file
        self.params['orn_volt_fn_base'] =  '%s/orn_volt_' % ( self.params['volt_folder'])
        self.params['mit_volt_fn_base'] =  '%s/mit_volt_' % ( self.params['volt_folder'])
        self.params['mit_glom_volt_fn_base'] =  '%s/mit_glom_volt_' % ( self.params['volt_folder'])
        self.params['mit_dend_volt_fn_base'] =  '%s/mit_dend_volt_' % ( self.params['volt_folder'])
        self.params['mit_prim_volt_fn_base'] =  '%s/mit_prim_volt_' % ( self.params['volt_folder'])
        self.params['gran_volt_fn_base'] =  '%s/gran_volt_' % ( self.params['volt_folder'])
        self.params['gran_periph_volt_fn_base'] =  '%s/gran_periph_volt_' % ( self.params['volt_folder'])
        self.params['gran_deep_volt_fn_base'] =  '%s/gran_deep_volt_' % ( self.params['volt_folder'])
        self.params['pg_volt_fn_base'] =  '%s/pg_volt_' % ( self.params['volt_folder'])
        self.params['pg_periph_volt_fn_base'] =  '%s/pg_periph_volt_' % ( self.params['volt_folder'])
        self.params['pg_deep_volt_fn_base'] =  '%s/pg_deep_volt_' % ( self.params['volt_folder'])
        self.params['pyr_volt_fn_base'] =  '%s/pyr_volt_' % ( self.params['volt_folder'])
        self.params['pyr_ca_fn_base'] =  '%s/pyr_calcium_' % ( self.params['volt_folder'])
        self.params['basket_volt_fn_base'] =  '%s/basket_volt_' % ( self.params['volt_folder'])
        self.params['rsnp_volt_fn_base'] =  '%s/rsnp_volt_' % ( self.params['volt_folder'])
        self.params['test_volt_fn_base'] =  '%s/test_volt_' % ( self.params['volt_folder'])
        self.params['test_isyn_fn_base'] =  '%s/test_isyn_' % ( self.params['volt_folder'])
        self.params['test_gsyn_fn_base'] =  '%s/test_gsyn_' % ( self.params['volt_folder'])
        self.params['readout_volt_fn_base'] =  '%s/readout_volt_spiking_' % ( self.params['volt_folder'])
        self.params['readout_mean_volt_fn'] =  '%s/readout_meanvolt.dat' % ( self.params['volt_folder'])
        self.params['readout_mean_volt_fig'] =  '%s/readout_meanvolt.png' % ( self.params['figure_folder'])

        # voltage files: each cell has its own file
        self.params['orn_iampa_fn_base'] =  '%s/orn_iampa_' % ( self.params['isyn_folder'])
        self.params['orn_igaba_fn_base'] =  '%s/orn_igaba_' % ( self.params['isyn_folder'])
        self.params['orn_iodor_fn_base'] =  '%s/orn_iodor_' % ( self.params['isyn_folder'])
        self.params['mit_iampa_fn_base'] =  '%s/mit_iampa_' % ( self.params['isyn_folder'])
        self.params['mit_igaba_fn_base'] =  '%s/mit_igaba_' % ( self.params['isyn_folder'])
        self.params['mit_inmda_fn_base'] =  '%s/mit_inmda_' % ( self.params['isyn_folder'])
        self.params['pg_iampa_fn_base'] =  '%s/pg_iampa_' % ( self.params['isyn_folder'])
        self.params['pg_igaba_fn_base'] =  '%s/pg_igaba_' % ( self.params['isyn_folder'])
        self.params['pg_inmda_fn_base'] =  '%s/pg_inmda_' % ( self.params['isyn_folder'])
        self.params['gran_iampa_fn_base'] =  '%s/gran_iampa_' % ( self.params['isyn_folder'])
        self.params['gran_igaba_fn_base'] =  '%s/gran_igaba_' % ( self.params['isyn_folder'])
        self.params['gran_inmda_fn_base'] =  '%s/gran_inmda_' % ( self.params['isyn_folder'])
        self.params['pyr_iampa_fn_base'] =  '%s/pyr_iampa_' % ( self.params['isyn_folder'])
        self.params['pyr_igaba_fn_base'] =  '%s/pyr_igaba_' % ( self.params['isyn_folder'])
        self.params['pyr_inmda_fn_base'] =  '%s/pyr_inmda_' % ( self.params['isyn_folder'])
        self.params['rsnp_iampa_fn_base'] =  '%s/rsnp_iampa_' % ( self.params['isyn_folder'])
        self.params['rsnp_igaba_fn_base'] =  '%s/rsnp_igaba_' % ( self.params['isyn_folder'])
        self.params['rsnp_inmda_fn_base'] =  '%s/rsnp_inmda_' % ( self.params['isyn_folder'])
        self.params['basket_iampa_fn_base'] =  '%s/basket_iampa_' % ( self.params['isyn_folder'])
        self.params['basket_igaba_fn_base'] =  '%s/basket_igaba_' % ( self.params['isyn_folder'])
        self.params['basket_inmda_fn_base'] =  '%s/basket_inmda_' % ( self.params['isyn_folder'])

        # nspike files: [nspikes, gid]
        # spiketimes files: [time, gid] (all cells in one file, but each process writes into a seperate file
        self.params['orn_spike_fn_base'] =  '%s/orn_nspikes_' % ( self.params['nspikes_folder'])
        self.params['orn_spikes_merged_fn_base'] =  '%s/orn_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['orn_spiketimes_fn_base'] =  '%s/orn_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['orn_spiketimes_merged_fn_base'] =  '%s/orn_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['mit_spike_fn_base'] =  '%s/mit_nspikes_' % ( self.params['nspikes_folder'])
        self.params['mit_spikes_merged_fn_base'] =  '%s/mit_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['mit_spiketimes_fn_base'] =  '%s/mit_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['mit_spiketimes_merged_fn_base'] =  '%s/mit_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['mit_response'] = '%s/mit_response_not_normalized_np%d' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_response_normalized'] = '%s/mit_response_normalized_np%d' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_rescaled'] = '%s/mit_nspikes_rescaled_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_cells'] = '%s/mit_nspikes_normed_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_patterns'] = '%s/mit_nspikes_normed_patterns_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_glom_cells'] = '%s/mit_nspikes_normed_glom_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_patterns_then_cells'] = '%s/mit_nspikes_normed_patterns_then_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_cells_then_patterns'] = '%s/mit_nspikes_normed_cells_then_patterns_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        # decide which mit response should be used as MDS input
#        self.params['mit_mds_input_fn'] = self.params['mit_response_normalized']
        self.params['mit_mds_input_fn'] = self.params['mit_nspikes_rescaled']

        self.params['gran_spike_fn_base'] =  '%s/gran_nspikes_' % ( self.params['nspikes_folder'])
        self.params['gran_spikes_merged_fn_base'] =  '%s/gran_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['gran_spiketimes_fn_base'] =  '%s/gran_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['gran_spiketimes_merged_fn_base'] =  '%s/gran_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['pg_spike_fn_base'] =  '%s/pg_nspikes_' % ( self.params['nspikes_folder'])
        self.params['pg_spikes_merged_fn_base'] =  '%s/pg_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['pg_spiketimes_fn_base'] =  '%s/pg_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['pg_spiketimes_merged_fn_base'] = '%s/pg_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['ob_spikes_merged_fn_base'] =  '%s/ob_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['ob_spiketimes_merged_fn_base'] = '%s/ob_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['pyr_spike_fn_base'] =  '%s/pyr_nspikes_' % ( self.params['nspikes_folder'])
        self.params['pyr_spikes_merged_fn_base'] =  '%s/pyr_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['pyr_spiketimes_fn_base'] =  '%s/pyr_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['pyr_spiketimes_merged_fn_base'] =  '%s/pyr_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['basket_spike_fn_base'] =  '%s/basket_nspikes_' % ( self.params['nspikes_folder'])
        self.params['basket_spikes_merged_fn_base'] =  '%s/basket_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['basket_spiketimes_fn_base'] =  '%s/basket_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['basket_spiketimes_merged_fn_base'] =  '%s/basket_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['rsnp_spike_fn_base'] =  '%s/rsnp_nspikes_' % ( self.params['nspikes_folder'])
        self.params['rsnp_spikes_merged_fn_base'] =  '%s/rsnp_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['rsnp_spiketimes_fn_base'] =  '%s/rsnp_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['rsnp_spiketimes_merged_fn_base'] =  '%s/rsnp_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['oc_spikes_merged_fn_base'] =  '%s/oc_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['oc_spiketimes_merged_fn_base'] = '%s/oc_spiketimes_merged_' % ( self.params['spiketimes_folder'])
        self.params['readout_spike_fn_base'] =  '%s/readout_nspikes_' % ( self.params['nspikes_folder'])
        self.params['readout_spikes_merged_fn_base'] =  '%s/readout_nspikes_merged_' % ( self.params['nspikes_folder'])
        self.params['readout_spiketimes_fn_base'] =  '%s/readout_spiketimes_' % ( self.params['spiketimes_folder'])
        self.params['readout_spiketimes_merged_fn_base'] =  '%s/readout_spiketimes_merged_' % ( self.params['spiketimes_folder'])
#        self.params['ob_output_poisson_fn_base'] ='%s/ob_output_poisson_' % (self.params[''])



    def check_folders(self):
        """
        Returns True if all folders exist, False otherwise
        """
        all_folders_exist = True
        for f in self.params['folder_names']:
            if not os.path.exists(f):
                all_folders_exist = False

        return all_folders_exist

    def create_folders(self):
        """
        Must be called from 'outside' this class before the simulation
        """

        for f in self.params['folder_names']:
            if not os.path.exists(f):
                print 'Creating folder:\t%s' % f
                os.system("mkdir -p %s" % (f))

    def load_params(self):
        """
        return the simulation parameters in a dictionary
        """
        return self.params


    def update_values(self, kwargs):
        for key, value in kwargs.iteritems():
            self.params[key] = value
        # update the dependent parameters
        self.set_filenames()


    def write_parameters_to_file(self, fn=None):
        if not (os.path.isdir(self.params['folder_name'])):
            print 'Creating folder:\n\t%s' % self.params['folder_name']
            self.create_folders()

        if fn == None:
            fn = self.params['params_fn_json']
        print 'Writing parameters to: %s' % (fn)
        output_file = file(self.params['params_fn_json'], 'w')
        d = json.dump(self.params, output_file)


    def hoc_export(self):
        '''
        Write all parameters to a file, executable by NEURON.
        - numeric values are simple written to file
        - strings need a strdef statement before and sprint statement.
            A regular sprint statement looks like: sprint(var_name, 'something_important')
        '''
        hoc_file_fn = self.params['hoc_file']
        lines = ''
        for p in self.params.keys():
            val = self.params.get(p)
#            print val, type(val)
            if (type(val) == type(1.0)):
                lines += '%s = %f\n' % (p, val)
            elif (type(val) == type(0)):
                lines += '%s = %d\n' % (p, val)
            elif (type(val) == type('string')):
                lines += 'strdef %s\n' % p
                lines += 'sprint(%s, \"%s\")\n' % (p, val)
#            elif (type(val) == type([])):
#                lines += str(val) # will not work because it's a stupid HOC file!
            else:
                print 'type(val)', type(val)
#                print 'Could not write \'%s\' or %s to file %s' % (p, val, hoc_file_fn)
                print 'Could not write \'%s\' to hoc file, as hoc does not understand: %s' % (p, type(val))

        if not os.path.exists(hoc_file_fn):
            os.system('touch %s' % hoc_file_fn)
#            print 'Hoc file %s does not exist!' % (hoc_file_fn)
        hoc_file = file(hoc_file_fn, 'w')
        hoc_file.write(lines)
        hoc_file.close()

    def print_cell_gids(self):
        print "ORNs: %d\t\t%d -\t%d\n" % (self.params['n_orn'], self.params['orn_offset'], self.params['orn_offset'] + self.params['n_orn'] - 1)
        print "MIT:  %d\t\t%d -\t%d\n" % (self.params['n_mit'], self.params['mit_offset'], self.params['mit_offset'] + self.params['n_mit'] - 1)
        print "GRAN: %d\t\t%d -\t%d\n" % (self.params['n_gran'], self.params['gran_offset'], self.params['gran_offset'] + self.params['n_gran'] - 1)
        print "PG:   %d\t\t%d -\t%d\n" % (self.params['n_pg'], self.params['pg_offset'], self.params['pg_offset'] + self.params['n_pg'] - 1)
        print "PYR:  %d\t\t%d -\t%d\n" % (self.params['n_pyr'], self.params['pyr_offset'], self.params['pyr_offset'] + self.params['n_pyr'] - 1)
        print "BASKET:%d\t\t%d -\t%d\n" % (self.params['n_basket'], self.params['basket_offset'], self.params['basket_offset'] + self.params['n_basket'] - 1)
        print "RSNP: %d\t\t%d -\t%d\n" % (self.params['n_rsnp'], self.params['rsnp_offset'], self.params['rsnp_offset'] + self.params['n_rsnp'] - 1)
        print "READOUT:%d\t%d -\t%d\n" % (self.params['n_readout'], self.params['readout_offset'], self.params['readout_offset'] + self.params['n_readout'] - 1)

