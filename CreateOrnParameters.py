import numpy as np
import numpy.random as rnd
import matplotlib.mlab as mlab

class CreateOrnParameters(object):

    def __init__(self, param_dict):

        self.params = param_dict
        self.orn_params_fn_base = self.params['orn_params_fn_base']
        self.n_or = self.params['n_or']
        self.n_gor = self.params['n_gor']
        self.n_orn_x = self.params['n_orn_x']
        self.n_patterns = self.params['n_patterns']
        self.n_orn_y = self.params['n_orn_y']
        self.n_orn = self.params['n_orn']
        self.num_params = 9 + 1 # + 1 for the id

        self.param_list = []
        self.min_values = []
        self.max_values = []

        self.gor_values = np.zeros(self.n_orn_x).tolist()

        self.activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))
        self.vgor = np.zeros(self.n_orn)
        self.voor = np.zeros(self.n_orn)
        self.vc_Kd= np.zeros(self.n_orn)
        self.vgna = np.zeros(self.n_orn)
        self.vgk = np.zeros(self.n_orn)
        self.vgkcag = np.zeros(self.n_orn)
        self.vgcal = np.zeros(self.n_orn)
        self.vgleak= np.zeros(self.n_orn)
        self.vtau_cadec = np.zeros(self.n_orn) # time constant for calcium removal

        self.set_gor_values()
        self.current_param_values = []
        # a list of lists containing for each row the following parameters:
                #"gna","gk","gkcag","gcal","gleak_orn", "tau_cadec"]


    def create_parameters_for_concentration_sweep(self, noisy_patterns=False):
        """
        orthogonal patterns means a single odorant is presented to the system
        """
        if noisy_patterns:
            # This requires that there already exists a odorant-OR-affinity matrix
            self.add_pattern_noise(self.params['OR_affinity_noise'])
        else:
            self.create_activation_matrix_for_concentration_sweep()

        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = self.activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
#                if (activation > 1.0):
#                    activation = 1.0
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1
        # stop new


    def create_single_odorant_patterns(self, noisy_patterns=False):
        """
        orthogonal patterns means a single odorant is presented to the system
        """
#        self.create_rnd_activation_matrix()
#        self.create_trinormal_affinities()
        if noisy_patterns: # This requires that there already exists a odorant-OR-affinity matrix
            self.add_pattern_noise(self.params['OR_affinity_noise'])
        else:
            self.create_single_odorant_activation_matrix()

        for pn in xrange(self.params["n_patterns"]):
            # calculate the number of receptors activated by the current odorant or pattern
            for OR in xrange(self.params['n_or']): # equal to a row in the orn array
                # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs)
                activation = self.activation_matrix[pn, OR]# * self.params['scale_factor_activation_matrix']
                for x in xrange(self.params['n_orn_x']):
                    orn_id = OR * self.params['n_orn_x'] + x
                    self.voor[orn_id] = activation #* conc

            self.set_gor_values()
            self.set_gna_values()
            self.set_gk_values()
            self.set_gkcag_values()
            self.set_gcal_values()
            self.set_gleak_values()
            self.set_tau_cadec_values()
            output_fn = self.orn_params_fn_base + "%d.dat" % (pn)
            print "Writin orn parameters for pattern \t... %d / %d to file: %s" % (pn, self.n_patterns, output_fn)
            self.write_params_to_file(output_fn)
        return 1


    def create_single_odorant_activation_matrix(self):
        """
        create_single_odorant_activation_matrix 
        """
        self.activation_matrix = np.zeros((self.params['n_patterns'], self.params['n_or']))

        """
        Distances are drawn from a tri-modal normal distributions
        This is done by choosing one of the three possible distributions 
        according to its probability p_i. 
        The probability with which one of the normal distributions is chosen 
        is determined by the respective integral (with respect to the other normal distributions).
        """

        p = self.params['odorant_receptor_distance_distribution_parameters']
        dist_range = self.params['odorant_receptor_distance_range']
        x = np.linspace(dist_range[0], dist_range[1], 1000)
        A1 = np.array(p[0] * mlab.normpdf(x, p[1], p[2])).sum() # integral of Gauss 1
        A2 = np.array(p[3] * mlab.normpdf(x, p[4], p[5])).sum() # integral of Gauss 2
        A3 = np.array(p[6] * mlab.normpdf(x, p[7], p[8])).sum() # integral of Gauss 3
        p1 = A1 / (A1 + A2 + A3)
        p2 = A2 / (A1 + A2 + A3)
        p3 = A3 / (A1 + A2 + A3)

        eps = 1e-12
        rnd.seed(self.params['seed_activation_matrix'])
        distances = np.zeros((self.params['n_patterns'], self.params['n_or']))
        for pn in xrange(self.params['n_patterns']):
            for OR in xrange(self.params['n_or']):
                # draw a random distance from the fitted distance distribution
                dist = self.odorant_odor_distance_distribution((p1, p2, p3))
                distances[pn, OR] = dist
                affinity = np.exp(-dist * self.params['distance_affinity_transformation_parameter'])
#                print 'DEBUG pn %d OR %d \tdist = %.6f\taffinity = %.6f' % (pn, OR, dist, affinity)
                self.activation_matrix[pn, OR] = affinity

        print 'Distances.mean:', distances.mean()
        print 'Distances.median:', np.median(distances)

        if self.params['OR_activation_normalization']:
#            for pn in xrange(self.params['n_patterns']):
#                self.activation_matrix[pn, :] /= self.activation_matrix[pn, :].max()
#                self.activation_matrix[pn, :] /= self.activation_matrix[pn, :].sum()

            for OR in xrange(self.params['n_or']):
                self.activation_matrix[:, OR] /= self.activation_matrix[:, OR].sum()

        print "Activation matrix fn:", self.params['activation_matrix_fn']
        np.savetxt(self.params['activation_matrix_fn'], self.activation_matrix)



    def odorant_odor_distance_distribution(self, which_gauss):
        """
        Returns a distance between a virtual OR and a virtual odorant.
        This can in principle be any distribution.
        Here, we use a tri-modal Gaussian distribution which has been fitted to
        the real distance distribution gained by k-means clustering (the centroids being the ORs).
        For details see script: average_OR_affinity_distributions.py 
        Keyword arguments:
        which_gauss -- tuple of size three containing the probabilities with which each gauss distributions is picked
                        (p1, p2, p3) = which_gauss
        """
        (p1, p2, p3) = which_gauss
        p = self.params['odorant_receptor_distance_distribution_parameters']
        dist_range = self.params['odorant_receptor_distance_range']
        which_gauss = np.random.uniform(0, 1.)
        if which_gauss < p1:
            print 'p1 ', 
            return np.random.normal(p[1], p[2])
        elif (which_gauss < p2 + p1):
            print 'p2 ', 
            return np.random.normal(p[4], p[5])
        elif (which_gauss < p3 + p2 + p1):
            print 'p3 ', 
            return np.random.normal(p[7], p[8])





    def create_params_for_response_curve(self):
        """
        This function creates the parameters needed to
        measure a response curve or concentration sweep
        (one odorant with maximum affinity and different concentrations
        """
#        output_fn = self.params["orn_params_test"]
        output_fn = self.orn_params_fn_base + "0.dat"
        print "Writing orn parameters to " , output_fn
        self.set_oor_value_for_conc_sweep()
        self.set_gor_values()
        self.set_gna_values()
        self.set_gk_values()
        self.set_gkcag_values()
        self.set_gcal_values()
        self.set_gleak_values()
        self.set_tau_cadec_values()
        self.write_params_to_file(output_fn)
        return 1


    def set_oor_value_for_conc_sweep(self):
        # c / Kd values are equally distributed on a log10 scale
        # c_Kd = 10**exp_min .. 10**exp_max
        exp_min = -2.0
#        exp_min = -1.5
        exp_max = 2.0
        z = (exp_max - exp_min) / float(self.n_or)
        c_Kd = []
        c_Kd.append(10**exp_min)
        for i in xrange(self.n_orn_y - 1):
            c_Kd.append(c_Kd[-1] * 10**z)

        # all ORNs within one row share the same Kd and concentration input value (because same OR is expressed by the these ORNs) and concentration increases with the row
        for y in xrange(self.n_orn_y):
            for x in xrange(self.n_orn_x): # increasing gor
                orn_id = y * self.n_orn_x + x
                self.vc_Kd[orn_id] = c_Kd[y]
                self.voor[orn_id] = 1 - 1. / (1 + c_Kd[y])


    def set_gor_values(self):

        def gor_func(x, a1, a2, exp):
            """this function calculates the gor for  a given x """
            return a1 * x**exp + a2
        # --------- CREATE ORN PARAMETERS -------------
        # gor values are distributed between a min and max value
        # gor(x) = a * x^gor_exp + b
        # gor(x_min) = b = gor_min
        # gor(x_max) = gor_max
        # parameters for gor_func, i.e. how gor is distributed among the n_gor different ORNs expressing the same OR
        x_min = 0
        x_max = self.n_orn_x

        gor_exp = self.params["gor_exp"]
        b = self.params["gor_min"]
        a = (self.params["gor_max"] - self.params["gor_min"]) / ((x_max - 1)**self.params["gor_exp"] - x_min**self.params["gor_exp"])
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                gor_value = gor_func(j, a, b, gor_exp)
                self.vgor[orn_id] = gor_value
                self.gor_values[j] = gor_value
        return self.gor_values


    def set_gna_values(self):
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function gna(gor)
                self.vgna[orn_id] = self.params["gna"]


    def set_gk_values(self):
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function gk(gor)
                self.vgk[orn_id] = self.params["gk"]


    def set_gkcag_values(self):
        # before this function is called set_gor_values must be called
        gkcag_values = self.linear_transformation(self.gor_values, self.params['gkcag_params'][0], self.params['gkcag_params'][1])

        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgkcag[orn_id] = gkcag_values[j]

                # map j into the interval [gor_min, gor_max] which is the x value for the function gkcag(gor)
#                x = self.vgor[orn_id]
#                gkcag = p[0] * np.log10(x) + p[1]
#                self.vgkcag[orn_id] = gkcag


    def set_gcal_values(self):
        # before this function is called set_gor_values must be called
        gcal_values = self.linear_transformation(self.gor_values, self.params['gcal_params'][0], self.params['gcal_params'][1])
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgcal[orn_id] = gcal_values[j]
                # map j into the interval [gor_min, gor_max] which is the x value for the function gcal(gor)
#                x = self.vgor[orn_id]
#                gcal = p[0] * np.log10(x) + p[1]
#                self.vgcal[orn_id] = gcal


    def set_gleak_values(self):
        # before this function is called set_gor_values must be called
        gleak_values = self.linear_transformation(self.gor_values, self.params['gleak_params'][0], self.params['gleak_params'][1])


        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                self.vgleak[orn_id] = gleak_values[j]

                # map j into the interval [gor_min, gor_max] which is the x value for the function gleak(gor)
#                x = self.vgor[orn_id]
#                gleak = p[0] * np.exp(p[1] * x ** 2) + p[2]
#                self.vgleak[orn_id] = gleak


    def set_tau_cadec_values(self):
        # before this function is called set_gor_values must be called
        for i in xrange(self.n_orn_y):  # y-axis
            for j in xrange(self.n_orn_x):
                orn_id = i * self.n_orn_x + j
                # map j into the interval [gor_min, gor_max] which is the x value for the function tau_cadec(gor)
                self.vtau_cadec[orn_id] = self.params["tau_cadec"]


    def write_params_to_file(self, output_fn):
        """
        Write the ORN parameters into a NEURON readable file
        """
        print "writing orn params to ", output_fn
        orn_pf = file(output_fn, 'w')
        num_params = 9 + 1 # + 1 for the id
        first_line = "%d %d\n"  % (self.n_orn, num_params) # neuron readability
        orn_pf.write(first_line)

        # write these values to file
        for row in xrange(self.params["n_orn_y"]):
            for col  in xrange(self.params["n_orn_x"]):
                orn_id = col + row * self.params["n_orn_x"]
                line = "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % \
                        (orn_id, self.voor[orn_id], self.vc_Kd[orn_id], self.vgor[orn_id], \
                            self.vgna[orn_id], self.vgk[orn_id], self.vgkcag[orn_id], self.vgcal[orn_id],\
                            self.vgleak[orn_id], self.vtau_cadec[orn_id])
                orn_pf.write(line)
        orn_pf.close()


    def write_current_param_values_to_file(self, output_fn=None):
        """
        This function is required for hand tuning the orn parameters.
        """
        if output_fn == None:
            output_fn = self.params['orn_params_fn_base'] + '0.dat'
        print "writing orn params to ", output_fn
        orn_pf = file(output_fn, 'w')
        num_params = 9 + 1 # + 1 for the id
        first_line = "%d %d\n"  % (self.n_orn, num_params) # neuron readability
        orn_pf.write(first_line)

        # write these values to file
        for row in xrange(self.params["n_orn_y"]):
            for col  in xrange(self.params["n_orn_x"]):
                orn_id = col + row * self.params["n_orn_x"]
                line = "%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n" % (orn_id, \
                        self.voor[orn_id], self.vc_Kd[orn_id], self.vgor[orn_id], \
                            self.current_param_values[col][0],\
                            self.current_param_values[col][1],\
                            self.current_param_values[col][2],\
                            self.current_param_values[col][3],\
                            self.current_param_values[col][4],\
#                            self.vgleak[orn_id],
                            self.current_param_values[col][5])
                orn_pf.write(line)
        orn_pf.flush()
        orn_pf.close()


    def linear_transformation(self, x, y_min, y_max):
        """
        x : the range to be transformed
        y_min, y_max : lower and upper boundaries for the range into which x
        is transformed to
        Returns y = f(x), f(x) = m * x + b
        """
        x_min = np.min(x)
        x_max = np.max(x)
        if x_min == x_max:
            x_max = x_min * 1.0001
        return (y_min + (y_max - y_min) / (x_max - x_min) * (x - x_min))


