import json
import numpy as np
import numpy.random as rnd
import os
import utils

class parameter_storage(object):
    """
    This class contains the simulation parameters in a dictionary called params.
    """

    def __init__(self):

        self.params = {}
        self.set_default_params()
        self.set_filenames()

    def set_default_params(self):
        pass


    def set_folder_name(self, folder_name=None):

        self.params['folder_name'] = 'Test'
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

        self.params['params_fn_json'] = '%ssimulation_parameters.json' % (self.params['parameters_folder'])
        self.params['hoc_file'] = '%s/simulation_params.hoc' % self.params['folder_name']
        self.params['info_file'] =  '%s/simulation_parameters_%d.info' % (self.params['folder_name'], self.params['sim_id']) # human readable file, format: parameter = value


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

