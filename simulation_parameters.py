import json
import numpy as np
import numpy.random as rnd
import os

class parameter_storage(object):
    """
    This class contains the simulation parameters in a dictionary called params.
    """

    def __init__(self, params_fn=None):
        """
        If a filename is given, it loads the json file and returns the dictionary.
        Else, it sets the default parameters.
        """

        self.params = {}
        if params_fn == None:
            self.set_default_params()
            self.set_filenames()
        else:
            self.load_params(params_fn)

        self.print_cell_gids()

    def set_default_params(self):


        self.params['OR_activation_normalization'] = True
        self.params['with_artificial_orns'] = 0
        
        self.params['Cluster'] = 1
        self.params['concentration_sweep'] = 0
        self.params['n_patterns'] = 50
        self.params['n_proc'] = 8   # on how many processors do you want to run the neuron code?
        self.params['ob_oc_random_conns'] = False
        self.params['oc_oc_random_conns'] = False

        self.params['with_noise'] = 1
        self.params['with_bias'] = 1 #if 0: you should remove insert gk_ka from the pyr_rs template! and recompile the neuron_files
        self.params['with_cond_bias'] = 1 # if 0: there's a constant negative current into some pyr cells

        # ------ S E E D S -------------
        self.params['seed_activation_matrix'] = 312
        self.params['seed'] = 9 # this is for pattern generation, weight randomization, etc
        self.params['seed_connections'] = 0
        self.params['netstim_seed'] = 10 # netstim_seed acts as an offset for the RNG provided to the noise input via NetStims

        # ------ N E T W O R K    S I Z E -----------
        if (self.params['concentration_sweep'] == 1):  # it's a concentration sweep
            self.params['n_patterns'] = 1
            # if concentration_sweep: n_or represents the number of different concentrations measured
            if (self.params['Cluster'] == 1):
                self.params['n_or'] = 32
            else:
                self.params['n_or'] = 24
        else:
#            self.params['n_or'] = 40
#            self.params['n_or'] = 40
            self.params['n_or'] = self.params['n_patterns']
        if (self.params['Cluster'] == 1):
            self.params['rel_orn_mit'] = 200
            self.params['rel_gran_mit'] = 100# number of granule cells per mitral cell
            self.params['rel_pg_mit']  = 20# number of periglomerular cells per mitral cell, ~ 20 according to Shepherd
        else:
            self.params['rel_orn_mit'] = 20
            self.params['rel_gran_mit'] = 10# number of granule cells per mitral cell
            self.params['rel_pg_mit']  = 10# number of periglomerular cells per mitral cell, ~ 20 according to Shepherd

        self.params['print_debug'] = 1 # flag to print more or less output
        # ------ C E L L     N U M B E R S ------------
#        self.params['n_gor'] = 16# number of mitral cells per glomerulus
        self.params['n_gor'] = 8# number of mitral cells per glomerulus
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
        self.params['n_test_mit'] = 10
        self.params['n_test_gran'] = 10
        self.params['n_test_pg'] = 20
        self.params['n_test_pyr'] = 5
        self.params['n_sample_pyr_per_mc'] = 1
        self.params['n_test_basket'] = 1
        self.params['n_sample_basket_per_hc'] = 1
#        self.params['n_sample_basket_per_hc'] = int(round(self.params['n_basket_per_hc'] / 10.0))
        self.params['n_test_rsnp'] = 1
        self.params['n_sample_rsnp_per_mc'] = 1

        # BCPNN parameters
        self.params['p_ij_thresh'] = 1e-8



        # ---------------- E X P E R I M E N T A L    P A R A M E T E R S --------- #
        self.params['temperature'] = 36# [Celsius] very important! required for NEURON simulations
        self.params['t_sim']	= 1600 # [ms] simulated time
        self.params['time_step']= 0.025   # [ms] max time step
        self.params['time_step_rec']= 0.5  # [ms] time step for recording membrane potentials etc
        self.params['thresh']	= 0     # [mV] threshold for spike detection. thresh is currently the same for all cells 	
        self.params['t_start']	= 0     # [ms] start time of current injection into orn cells
        self.params['tau_odorinput_sigmoid'] = 20 # [ms] time constant for sigmoidal function for input conductance time course, check with neuron_files/odorinput.mod
        self.params['t_stop']	= 25 * self.params['tau_odorinput_sigmoid'] # [ms] start time for decaying sigmoid for odor input conductance
#        self.params['curr_amp']= 100	# [nA] amplitude of current injected into orn cells
        self.params['v_init'] = -70.	# [mV]
        self.params['v_init_sigma'] = 5 # [mV]



        # ODORANT - OR DISTRIBUTION PARAMERS
        # obtained through average_OR_affinity_distributions.py
        # The odorant_receptor_distance_range marks the range of possible distances between ORs and odorants based
        # on the clustering results obtained from average_OR_affinity_distributions.py
        self.params['odorant_receptor_distance_range'] = (0, 4.330310991999920844e+01)
        self.params['odorant_receptor_distance_distribution_parameters'] = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
         1.909487e+01, 1.110809e+01, 3.353855e+00, \
         4.188897e+00, 4.088460e+01, 4.966478e-01] # these values are taken from clustering the odorant space with 40 ORs

        # ---------------- C E L L    P A R A M E T E R S --------- # 
        # ---------------- ORN cell parameters:

        # gor stands for the maximum conductance evoked by an odor
        # gor values are distributed between a min and max value
        # good values for system without noise
#        self.params['gor_min'] = 3e-5 
#        self.params['gor_max'] = 5e-4
        # if ORNs have all the same conductance parameters, this is the list:
        self.params['gna'] = 0.5        # [S/cm2]
        self.params['gk'] = 0.05        # [S/cm2]
        self.params['gkcag'] = 0.01     # [S/cm2]
        self.params['gcal'] = 6e-4      # [S/cm2]
        self.params['gleak_orn'] = 1e-4 # [S/cm2]
        self.params['tau_cadec'] = 1000 # [ms]

        # parameters for gleak, gkcag, gcal, gained through combined hand tuning / fitting procedure 
        self.params['gor_params'] = [5e-5, 1.1e-3]
        self.params['gor_min'] = self.params['gor_params'][0]
        self.params['gor_max'] = self.params['gor_params'][1]
        self.params['gor_exp'] = 2
        self.params['gkcag_params'] = [5e-3, 5e-2]
        self.params['gcal_params'] =  [3e-5, 0.8e-5]
        self.params['gleak_params'] = [8.0e-5, 1.2e-4]
#        self.params['gkcag_params'] = [5e-3, 5e-2]
#        self.params['gcal_params'] =  [1e-5, 1e-5]
#        self.params['gleak_params'] = [1.2e-4, 8e-5]

#        self.params['gkcag_params'] = [0.005000, 0.005000]
#        self.params['gcal_params'] =  [0.000050, 0.000500]
#        self.params['gleak_params'] = [0.000050, 0.000300]

#        self.params['gkcag_params'] = [4.99086530e-03, 2.26738160e-02, 2.26738160e-02]
#        self.params['gcal_params'] =  [4.99086531e-04, 2.26738160e-03]
#        self.params['gleak_params'] = [4.25453912e-05, -5.18713818e+05, 3.47077557e-05]

#        self.params['gkcag_params'] = [4.99086530e-03, 2.26738160e-02, 2.26738160e-02]
#        self.params['gcal_params'] =  [4.99086531e-04, 2.26738160e-03, 2.26738160e-03]
#        self.params['gcal_params'] =  [4.99086531e-04, 2.0e-3]
#        self.params['gleak_params'] = [5.0e-05, -6.0e+05, 1.0e-05]


        # --------------- Artificial ORNs
        self.params['f_out_artificial_orn_min'] = 30
        self.params['f_out_artificial_orn_max'] = 40
        self.params['t_sigma_artificial_orn'] = 50 # [ms] -> what's the temporal spread of spikes originating from ORNs?


        # ---------------- OB connectivity parameters
        # ---------------- ORN -> MIT connectivity
        self.params['w_nmda_mult'] = 3 # ORN - MT and ORN - PG: NMDA weights are multiplied by this factor compared to AMPA weights
#        self.params['w_nmda_mult'] = 2 # ORN - MT and ORN - PG: NMDA weights are multiplied by this factor compared to AMPA weights
        self.params['with_auto_receptors'] = 1 # flag for glutamatergic autoreceptors on mitral cells
        self.params['w_mit_ampa_autoreceptors'] = 0.005 # weight of the NetCons in the mitral cell primary dendrite representing AMPA autoreceptors
#        self.params['w_mit_ampa_autoreceptors'] = 0.002 # weight of the NetCons in the mitral cell primary dendrite representing AMPA autoreceptors
        self.params['w_mit_nmda_autoreceptors'] = self.params['w_mit_ampa_autoreceptors'] * self.params['w_nmda_mult'] # weight of the NetCons in the mitral cell primary dendrite representing NMDA autoreceptors
        self.params['w_orn_mit_target'] = 0.1 # target excitatory conductance received by a mitral cell
        self.params['w_orn_mit_sigma'] = 0.1 # sigma of the normal distribution for drawing conn weights
        self.params['w_orn_mit_mult'] = 6.0 # orns with lower sensitivity have smaller output rates at high concentrations, thus their outgoing connection weight to MT and PG cells is multiplied by this factor

#        self.params['w_orn_mit_mult'] = 1.0 # orns with lower sensitivity have smaller output rates at high concentrations, thus their outgoing connection weight to MT and PG cells is multiplied by this factor

        # the weight from some ORN groups to their target MIT is multiplied to compensate for their lower output rate
        self.params['orn_mit_change_ids'] = range(8)
        self.params['orn_mit_change_factors'] = [1.25, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#        self.params['orn_mit_change_factors'] = [2.0, 2.0, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0]
#        self.params['orn_mit_change_factors'] = [1.5, 1.5, 1.4, 1.3, 1.1, 1.0, 1.0, 1.0]
#        self.params['orn_mit_change_factors'] = [2.1, 1.9, 1.6, 1.5, 1.2, 1.1, 1.0, 1.0]
#        self.params['orn_mit_change_ids'] = [0, 1, 2, 3, 4]
#        self.params['orn_mit_change_factors'] = [1.9, 1.7, 1.4, 1.3, 1.1]
        assert (len(self.params['orn_mit_change_ids']) == len(self.params['orn_mit_change_factors']))
        self.params['orn_pg_change_ids'] = range(8)# index of PG cells whose orn-pg weights are not increased (=index of a mitral cell with too low f_out in the interval code response curve plot
#        self.params['orn_pg_change_factors'] = [2.0, 2.0, 2.0, 1.5, 1.5, 1.5, 1.0, 1.0]
#        self.params['orn_pg_change_factors'] = [2.0, 2.0, 2.0, 1.5, 1.5, 1.5, 1.0, 1.0]
        self.params['orn_pg_change_factors'] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
#        self.params['orn_pg_change_factors'] = [8.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.0, 1.0]
        assert (len(self.params['orn_pg_change_factors']) == len(self.params['orn_pg_change_ids']))

        # ---------------- ORN -> PG connectivity
        self.params['w_orn_pg_target'] = 0.008 # Sum of excitatory weights to be received by a PG cell
        self.params['w_orn_pg_sigma'] = 0.1 # sigma of the normal distribution for drawing conn weights
        self.params['w_orn_pg_mult'] = self.params['w_orn_mit_mult'] # orns with lower sensitivity have smaller output rates at max concentration, thus their outgoing connection weight to MT and PG cells is multiplied by this factor
        self.params['orn_inh_shift'] = 1.0 # inhibitory gaussian curve is shifted to the left by this number of n_orn_exc_sigma

        # MT o--- PG: serial
        self.params['w_pg_mit_serial_target'] = 8.0 # total inhibitory conductance to be received by one MT cell from the population of 'serial' periglomerular cells
#        self.params['w_pg_mit_serial_target'] = 5.0 # total inhibitory conductance to be received by one MT cell from the population of 'serial' periglomerular cells
        self.params['w_pg_mit_serial_sigma'] = 0.1
        # MT o---< PG: reciprocal
#        self.params['w_pg_mit_reciprocal_target'] = 0.004 # total inhibitory conductance to be received by one MT cell from the population of 'reciprocal' periglomerular cells
        self.params['w_pg_mit_reciprocal_target'] = 0.2 # total inhibitory conductance to be received by one MT cell from the population of 'reciprocal' periglomerular cells
        self.params['w_pg_mit_reciprocal_sigma'] = 0.1
        self.params['w_mit_pg_serial'] = 0.006 # total excitatory conductance to be received by one PG cell from one MT cell
        self.params['w_mit_pg_serial_sigma'] = 0.1
        self.params['w_mit_pg_reciprocal'] = 1e-3 # total excitatory conductance to be received by one PG cell from one MT cell
        self.params['w_mit_pg_reciprocal_sigma'] = 0.1
        # MT o----< GRAN: local
        self.params['n_mit_gran_syn_local'] = 2000 # number of reciprocal synapses between one MT cell and all granule cells within the same glomerulus
        self.params['w_mit_gran_local_target'] = 0.03 # total excitation received in average by a granule cell through local excitatory synapses from all MT cells within this glomerulus, ampa weights
#        self.params['w_mit_gran_local_target'] = 0.015 # total excitation received in average by a granule cell through local excitatory synapses from all MT cells within this glomerulus, ampa weights
#        self.params['w_mit_gran_local_target'] = 0.02 # total excitation received on average by a granule cell through local excitatory synapses from all MT cells within this glomerulus, ampa weights
        self.params['w_mit_gran_nmda_mult'] = 3 # dendro-dendritic synapses from MT onto Gran cells are dominantly NMDA mediated (Schoppa'98), thus multiply their weight compared to AMPA weights
        self.params['w_mit_gran_local_sigma'] = 0.1   # std deviation for gaussian distributed weights
        self.params['w_gran_mit_local_target'] = 10.0
#        self.params['w_gran_mit_local_target'] = 2.0
        self.params['w_gran_mit_local_sigma'] = 0.1
        # MT o----< GRAN: global
        if (self.params['concentration_sweep'] == 1):
            self.params['n_mit_gran_syn_global'] = 1 # number of reciprocal synapses between one MT cell and granule cells in other glomeruli
            self.params['w_mit_gran_global_target'] = 1e-9  # total excitatory conductance received by a Gran cells from non-local MT cells, i.e. MT cells in other glomeruli
            self.params['w_gran_mit_global_target'] = 1e-9  # total inhibitory conductance received by an MT cells from non-local DDI connections with Gran cells, i.e. Gran cells in other glomeruli
        else:
            self.params['n_mit_gran_syn_global'] = 100 # number of reciprocal synapses between one MT cell and granule cells in other glomeruli
            self.params['w_mit_gran_global_target'] = 0.002  # total excitatory conductance received by a Gran cells from non-local MT cells, i.e. MT cells in other glomeruli
            self.params['w_gran_mit_global_target'] = 0.1 # total inhibitory conductance received by an MT cells from non-local DDI connections with Gran cells, i.e. Gran cells in other glomeruli
        self.params['w_mit_gran_global_sigma'] = 0.1
        self.params['w_gran_mit_global_sigma'] = 0.1
        # --------------- OB parameters for dendro-dendritic inhibition
        self.params['ddi_thresh'] = -40
        self.params['ddi_mit_glom_thresh'] = self.params['ddi_thresh']
        self.params['ddi_mit_dend_thresh'] = self.params['ddi_thresh']
        self.params['autorec_mit_prim_thresh'] = self.params['ddi_thresh']
        self.params['autorec_mit_dend_thresh'] = self.params['ddi_thresh']
        self.params['ddi_pg_periph_thresh'] = self.params['ddi_thresh']
        self.params['ddi_gran_periph_thresh'] = self.params['ddi_thresh']
        self.params['ddi_pg_mit_delay'] = 1
        self.params['ddi_mit_pg_delay'] = 1
        self.params['ddi_mit_gran_delay'] = 1
        self.params['ddi_gran_mit_delay'] = 1
        self.params['mit_autoreceptor_delay'] = 1


        # --------------- CORTICAL CONNECTIVITY
        # all weights are given in uS, thus w=0.001 is 1 nS
        # within one minicolumn # [E_psp_height in mV at V_rest for pyr_rs]
        # a weight of 0.005 = EPSP_height = 4.25 mV
        self.params['w_pyr_pyr_local'] = 0.001
        self.params['w_pyr_basket'] = 0.003
        self.params['w_basket_pyr'] = 0.005
        self.params['w_basket_basket'] = 0.002
        self.params['w_rsnp_pyr'] = 0.002       # -0.8 mV
        self.params['w_nmda_mult_oc'] = 1.0     # for oc - oc connections
        self.params['w_pyr_readout'] = 0.001

        # pyr->pyr global:
#        self.params['w_pyr_pyr_global_max'] = 1e-8
#        self.params['w_pyr_rsnp_max'] = 1e-8
        self.params['w_pyr_pyr_global_max'] = 0.0005
        self.params['w_pyr_rsnp_max'] = 0.0005

        # weight variation
        self.params['w_sigma'] = 0.1 # 0.2
        self.params['w_pyr_pyr_local_sigma'] = self.params['w_sigma'] * self.params['w_pyr_pyr_local']
        self.params['w_pyr_pyr_global_sigma'] =self.params['w_sigma'] * self.params['w_pyr_pyr_global_max']
        self.params['w_pyr_basket_sigma'] =self.params['w_sigma'] * self.params['w_pyr_basket']
        self.params['w_pyr_rsnp_sigma'] = self.params['w_sigma'] * self.params['w_pyr_rsnp_max']
        self.params['w_basket_pyr_sigma'] = self.params['w_sigma'] * self.params['w_basket_pyr']
        self.params['w_basket_basket_sigma'] = self.params['w_sigma'] * self.params['w_basket_basket']
        self.params['w_rsnp_pyr_sigma'] = self.params['w_sigma'] * self.params['w_rsnp_pyr']

        # weight threshold: when drawing connections only weight/ bigger than this are finally drawn
        self.params['weight_threshold'] = 1e-5

        # connection probabilities
        self.params['p_rsnp_pyr'] = 0.7
        self.params['p_pyr_pyr_local'] = 0.25
        self.params['p_pyr_pyr_global'] = 0.3
        self.params['p_pyr_basket'] = 0.7
        self.params['p_basket_pyr'] = 0.7
        self.params['p_pyr_rsnp'] = 0.3
        self.params['p_basket_basket'] = 0.7

        # ---------------- MIT -> PYR connectivity
        self.params['w_mit_pyr_max'] = 0.005             # max weight for exc mit -> pyr connection
        self.params['w_mit_rsnp_max'] = 0.002             # max weight for exc mit -> rsnp connection, new 

#        self.params['w_ampa_thresh'] = 0.002            # weights (transformed to the detailed model) larger than this value will be connected also via an AMPA
        self.params['w_mit_pyr_sigma_frac'] = self.params['w_sigma']
        self.params['w_mit_rsnp_sigma_frac'] = self.params['w_mit_pyr_sigma_frac']


        # ---------------- CORTICAL CELLS:
        # Potassium M-current (adaptation strength)
        self.params['g_m_pyr'] = 8e-5 # [S / cm2]
        self.params['g_m_basket'] = 3e-5 # [S / cm2]
        self.params['g_m_rsnp'] = 4e-5 # [S / cm2]
        self.params['g_leak_readout'] = 8e-5 # [S / cm2]
        self.params['tau_max_g_m'] = 1000# [ms]
        self.params['g_ka_pyr_max'] = 50.0 # 40 # [uS / cm2] # bias conductance
        self.params['g_ka_readout_max'] = self.params['g_ka_pyr_max']
        self.params['i_bias_pyr_max'] = -0.15 # [pA] # pyramidal cells with the maximal bias value get this as negative iclamp.amp 
        self.params['i_bias_readout_max'] = -0.15 # [pA] # pyramidal cells with the maximal bias value get this as negative iclamp.amp 

        # Calcium gated potassium current
        self.params['g_kcag_pyr'] = 1e-5 # [S / cm2]
        # High threshold Calcium current
        self.params['g_cal_pyr'] = 1e-5 # [S / cm2]


        # ---------------- S Y N A P S E     P A R A M E T E R S --------- # 
        self.params['tau_syn_exc'] = 10 # [ms] AMPA synapses onto somota of Pyr, Basket and Rsnp cells
        self.params['tau_syn_inh'] = 20 # [ms] GABBA ergic synapses onto somata of Pyr, Basket and Rsnp cells
        self.params['tau_nmda'] = 150 # [ms] this is not passed to the NMDA synapse itself, it's estimated from a fit of an exponential to the conductance time course and only used in the ConductanceCalculator class

        


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


        # -------- MDS - VQ - BCPNN  Parameters ---------------
        self.params['vq_ob_oc_overlap'] = 6 # if vq_overlap == 0: only one target Hypercolumn per mitral cell
        self.params['n_bcpnn_steps'] = 1
        self.params['vq_oc_readout_overlap'] = 0 # if vq_overlap == 0: only one target Hypercolumn per mitral cell



    def set_folder_name(self, folder_name=None):
        """
        This function is called from set_filenames in order to update all filenames 
        with the given folder_name
        Keyword arguments:
        folder_name -- string
        """

#        folder_name = 'FullSystemTest_0'
        folder_name = 'FullSystemTest_np50_postLearning'
#        folder_name = 'FullSystemTest_np%d_normalized_OrnAct' % (self.params['n_patterns'])
#        folder_name = 'ObTest93_seed%d%d'% (self.params['seed'], self.params['netstim_seed'])

#        folder_name = 'ResponseCurvesEpthOb_6'
        if self.params['Cluster']:
            folder_name = 'Cluster_' + folder_name
            use_abspath = False
        else:
            use_abspath = True
            
        if use_abspath:
            self.params['folder_name'] = os.path.abspath(folder_name)
        else:
            self.params['folder_name'] = folder_name
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
        self.params['conn_list_pg_mit_serial'] =  '%s/conn_list_pg_mit_serial.dat' % ( self.params['conn_folder'])
        self.params['conn_list_pg_mit_reciprocal'] =  '%s/conn_list_pg_mit_reciprocal.dat' % ( self.params['conn_folder']) # contains only PG ---o MT
        self.params['conn_list_mit_pg_reciprocal'] =  '%s/conn_list_mit_pg_reciprocal.dat' % ( self.params['conn_folder']) # contains MT ---< PG
        self.params['conn_list_mit_pg_serial'] =  '%s/conn_list_mit_pg_serial.dat' % ( self.params['conn_folder']) # contains MT ---< PG
        self.params['conn_list_mit_gran_local'] =  '%s/conn_list_mit_gran_local.dat' % ( self.params['conn_folder'])
        self.params['conn_list_mit_gran_global'] =  '%s/conn_list_mit_gran_global.dat' % ( self.params['conn_folder'])
        self.params['conn_list_gran_mit_local'] =  '%s/conn_list_gran_mit_local.dat' % ( self.params['conn_folder'])
        self.params['conn_list_gran_mit_global'] =  '%s/conn_list_gran_mit_global.dat' % ( self.params['conn_folder'])
        self.params['conn_list_mit_pyr'] =  '%s/conn_list_mit_pyr.dat' % ( self.params['conn_folder']) # generated from Anders' abstract weight matrices
        self.params['conn_list_mit_rsnp'] =  '%s/conn_list_mit_rsnp.dat' % ( self.params['conn_folder']) # generated from Anders' abstract weight matrices
        self.params['conn_list_mit_oc'] =  '%s/conn_list_mit_oc.dat' % ( self.params['conn_folder'])
        self.params['conn_list_mit_oc_base'] =  '%s/conn_list_mit_oc_' % (self.params['conn_folder'])
        self.params['conn_list_layer23'] =  '%s/conn_list_layer23.dat' % ( self.params['conn_folder']) # pyr->basket, basket->pyr, rsnp->pyr, pyr->pyr within one MC


        # the following two files are created by the GetConnectionsQuickFix class based on the learning output 
        self.params['conn_list_pyr_mit'] =  '%s/conn_list_pyr_mit.dat' % ( self.params['conn_folder']) # inhibitory feedback connections from OC to OB
        self.params['conn_list_pyr_pyr_base'] =  '%s/conn_list_pyr_pyr_' % ( self.params['conn_folder']) # when parallel connection creation method is used, files are written with this basename
        self.params['conn_list_pyr_pyr'] =  '%s/conn_list_pyr_pyr.dat' % ( self.params['conn_folder']) # additional recurrent connections in OC
        self.params['conn_list_pyr_rsnp_base'] =  '%s/conn_list_pyr_rsnp_' % ( self.params['conn_folder'])  # when parallel connection creation method is used, files are written with this basename
        self.params['conn_list_pyr_rsnp'] =  '%s/conn_list_pyr_rsnp.dat' % ( self.params['conn_folder']) # pyr->rsnp connections 
        self.params['conn_list_pyr_readout'] =  '%s/conn_list_pyr_readout.dat' % ( self.params['conn_folder']) # pyr -> readout layer



        # Files for MDS - VQ - BCPNN
        self.params['binary_oc_activation_fn'] = '%s/binary_oc_activation.dat' % (self.params['other_folder']) # this file is created by MDSVQ.create_mitral_response_space and stores the initial activation of minicolumns to begin the BCPNN learning
        self.params['silent_mit_fn'] = '%s/silent_mitral_cells.txt' % (self.params['other_folder'])
        # machine learning output files
        self.params['mds_ob_oc_output_fn'] = '%s/mds_ob_oc_output.dat' % (self.params['other_folder']) # this file stores the coordinates for the mitral cells in the mutual information space
        self.params['vq_ob_oc_output_fn'] = '%s/vq_ob_oc_output.dat' % (self.params['other_folder']) # this file stores the binary mit - hc connection matrix, (mit_hc_mask) created after VQ in the mutual information MDS space
        self.params['mit_mc_vq_distortion_fn'] = '%s/mit_mc_vq_distortion_' % (self.params['other_folder'])
        self.params['abstract_binary_conn_mat_ob_oc_fn'] = '%s/binary_conn_mat_ob_oc.dat' % (self.params['other_folder'])
        self.params['mds_oc_readout_output_fn'] = '%s/mds_oc_readout_output.dat' % (self.params['other_folder'])
        self.params['vq_oc_readout_output_fn'] = '%s/vq_oc_readout_output.dat' % (self.params['other_folder'])
        self.params['oc_readout_conn_fn'] = '%s/oc_readout_conn.dat' % (self.params['other_folder'])
        self.params['abstract_binary_conn_mat_oc_readout_fn'] = '%s/binary_conn_mat_oc_readout.dat' % (self.params['other_folder'])

        # cell_type to be attached after _fn_base
        self.params['mutual_information_fn_base'] = '%s/mutual_information_' % (self.params['other_folder'])# the actual mutual information between cells of a certain celltype, mi
        self.params['joint_entropy_fn_base'] = '%s/joint_entropy_' % (self.params['other_folder']) # the joint entropy between cells of a certain type, je
        self.params['information_distance_fn_base'] = '%s/information_distance_' % (self.params['other_folder']) # the information distance d = 1 - mi / je

        # spiking cortical activity, clustered by minicolumns
        self.params['clustered_mc_output'] = '%s/clustered_oc_output.dat' % self.params['figure_folder']
        self.params['clustered_mc_output_timebinned_fn_base'] = '%s/clustered_oc_output_timebinned_' % self.params['figure_folder']
        self.params['clustered_mc_output_wta'] = '%s/clustered_oc_output_wta.dat' % self.params['figure_folder']
        self.params['clustered_mc_output_normalized'] = '%s/clustered_oc_output_normed.dat' % (self.params['figure_folder'])

        # filenames for optional mds output of 2nd VQ in mitral cell response space
        self.params['mit_response_space_fn_base'] = '%s/mit_response_space' % self.params['other_folder'] # this file contains the 3D coordinates of mitral cells projecting to the same HC
        self.params['mit_response_space_centroids_fn_base'] = '%s/mit_response_space_centroids' % self.params['other_folder']# this file contains the 3D coordinates of mitral cells projecting to the same HC

        # filenames for BCPNN results
        # ob - oc
        self.params['oc_abstract_activity_fn'] = '%s/oc_abstract_activity.dat' % (self.params['bcpnn_folder'])
        self.params['ob_oc_abstract_weights_fn'] = '%s/ob_oc_abstract_weights.dat' % (self.params['bcpnn_folder'])
        self.params['ob_oc_abstract_bias_fn'] = '%s/ob_oc_abstract_bias.dat' % (self.params['bcpnn_folder'])
        # oc - oc : recurrent 
        self.params['oc_rec_abstract_activity_fn'] = '%s/oc_rec_abstract_activity.dat' % (self.params['bcpnn_folder'])
        self.params['oc_oc_abstract_weights_fn'] = '%s/oc_oc_abstract_weights.dat' % (self.params['bcpnn_folder'])
        self.params['oc_oc_abstract_bias_fn'] = '%s/oc_oc_abstract_bias.dat' % (self.params['bcpnn_folder'])
        # oc - readout
        self.params['readout_abstract_activity_fn']  = '%s/readout_abstract_activity.dat' % (self.params['bcpnn_folder'])
        self.params['oc_readout_abstract_weights_fn'] = '%s/oc_readout_abstract_weights.dat' % (self.params['bcpnn_folder'])
        self.params['oc_readout_abstract_bias_fn'] = '%s/oc_readout_abstract_bias.dat' % (self.params['bcpnn_folder'])



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
#            MIT - response - normalized:
#                a) the sum of spikes fired by each cell during all patterns is normalized to 1 -> each mitral cell has a normalized activty
#                b) if the sum of normalized activity of mitral cells within one glomerular unit > 1 -> set it to one
        self.params['mit_nspikes_rescaled'] = '%s/mit_nspikes_rescaled_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])  # Rescaled mit_spikes so that the global maximum = 1
        self.params['mit_nspikes_normed_cells'] = '%s/mit_nspikes_normed_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_patterns'] = '%s/mit_nspikes_normed_patterns_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_glom_cells'] = '%s/mit_nspikes_normed_glom_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_patterns_then_cells'] = '%s/mit_nspikes_normed_patterns_then_cells_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        self.params['mit_nspikes_normed_cells_then_patterns'] = '%s/mit_nspikes_normed_cells_then_patterns_np%d.dat' % (self.params['nspikes_folder'], self.params['n_patterns'])
        # decide which mit response should be used as MDS input
        self.params['mit_mds_input_fn'] = self.params['mit_response_normalized']
#        self.params['mit_mds_input_fn'] = self.params['mit_nspikes_rescaled'] 

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

        # connection matrices to be created after the training, abstract and detailed matrices for OB->OC, OC<->OC, OC->Readout
        self.params['connection_matrix_abstract_ob_oc_dat'] = '%s/connection_matrix_abstract_ob_oc.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_ob_oc_dat_unscaled'] = '%s/connection_matrix_abstract_ob_oc_unscaled.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_oc_dat'] = '%s/connection_matrix_abstract_oc_oc.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_oc_dat_unscaled'] = '%s/connection_matrix_abstract_oc_oc.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_ob_dat'] = '%s/connection_matrix_abstract_oc_ob.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_readout_dat'] = '%s/connection_matrix_abstract_oc_readout.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_ob_oc_unscaled_fig'] = '%s/connection_matrix_abstract_ob_oc_unscaled.png' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_ob_oc_fig'] = '%s/connection_matrix_abstract_ob_oc.png' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_oc_fig'] = '%s/connection_matrix_abstract_oc_oc.png' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_ob_fig'] = '%s/connection_matrix_abstract_oc_ob.png' % (self.params['conn_folder'])
        self.params['connection_matrix_abstract_oc_readout_fig'] = '%s/connection_matrix_abstract_oc_readout.png' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_ob_oc_dat'] = '%s/connection_matrix_detailed_ob_oc.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_oc_dat'] = '%s/connection_matrix_detailed_oc_oc.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_ob_dat'] = '%s/connection_matrix_detailed_oc_ob.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_readout_dat'] = '%s/connection_matrix_detailed_oc_readout.dat' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_ob_oc_fig'] = '%s/connection_matrix_detailed_ob_oc.png' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_oc_fig'] = '%s/connection_matrix_detailed_oc_oc.png' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_ob_fig'] = '%s/connection_matrix_detailed_oc_ob.png' % (self.params['conn_folder'])
        self.params['connection_matrix_detailed_oc_readout_fig'] = '%s/connection_matrix_detailed_oc_readout.png' % (self.params['conn_folder'])



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

    def load_params(self, fn):
        """
        Load a json-parameter file and 
        return the simulation parameters in a dictionary
        """
        f = file(fn, 'r')
        print 'Loading parameters from', fn
        self.params = json.load(f)
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
        d = json.dump(self.params, output_file, indent=0)


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

