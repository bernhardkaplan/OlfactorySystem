"""
This script analyses the results of cluster_odorant_space.py

cluster_odorant_space.py writes for each number of clusters and each trial one
distance matrix.

In order to compare, the average distribution of affinities (affinity-histogram, i.e.count vs. affinity value) are compared.
"""
import numpy as np
import os
import re
import pylab
import matplotlib.mlab as mlab
from scipy.optimize import leastsq
from scipy.special import erf # error function
import scipy.stats

# ---------- a bunch of distributions to try ------------------------------- 
def residuals_gauss(p, y, x):
    """
    x: normal x coordinate
    p: parameters of the function to fit, e.g. a and b in y = a * x + b
    """
    y1 = peval_gauss(x, p)
    err = y - y1
    return err 

def residuals_skew_normal(p, y, x):
    return y - peval_skew_normal(x, p)


def residuals_binormal(p, y, x):
    return y - peval_binormal(x, p)

def peval_skew_normal(x, p):
    # p[3] = amplitude
    # p[2] = alpha (skewness parameter)
    # p[1], p[0] = mu, sigma
    return p[3] * skew_normal(x, p[0], p[1], p[2])

def peval_gauss(x, p):
    return p[2] * mlab.normpdf(x, p[0], p[1])

def peval_binormal(x, p):
    # p[0] = w1
    # p[1] = mu1
    # p[2] = sigma1
    # p[3] = w2
    # p[4] = mu2
    # p[5] = sigma2
    return (p[0] * mlab.normpdf(x, p[1], p[2]) + p[3] * mlab.normpdf(x, p[4], p[5]))

def my_gauss(x, mu, sigma):
    return 1./(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5 * (1./sigma*(x - mu))**2)

def skew_normal(x, mu, sigma, alpha):
    # alpha  = skewness parameter
    return (1. / sigma) * mlab.normpdf((x - mu)/sigma, 0., 1.) * (1. + erf(alpha * (x - mu) / (sigma * np.sqrt(2))) )

def bimodal_normal(x, w1, mu1, sigma1, w2, mu2, sigma2):
    return w1 * mlab.normpdf(x, mu1, sigma1) + w2 * mlab.normpdf(x, m2, sigma2)
#    return 1./(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5 * (1./sigma*(x - mu))**2)


if __name__ == '__main__':

    folder = "OR_placement/AffinityMatrices/"
    figure_folder = 'OR_placement/Figures/'
    if not os.path.exists(figure_folder):
        os.system('mkdir %s' % figure_folder)
    folder = os.path.abspath(folder)

    n_bins = 100
    filter_fn = "affinity_matrix_OR"
    global_max = 0.
    # get the distribution of distances between centroids and virtual ORs
    # for each number of ORs
    for n_cluster in xrange(2, 20):
        file_cnt = 0
        # find the files matching the filter_fn pattern
        for fn in os.listdir(folder):
            m = re.match("%s%d_(\d+)\.dat" % (filter_fn, n_cluster), fn)
            if m:
                file_cnt += 1
        print "Averaging data for %d clusters, found %d files" % (n_cluster, file_cnt)

        # 1) get all the data written by cluster_odorant_space
        mean_count = np.zeros(n_bins)
        std_count = np.zeros(n_bins)
        all_data = np.zeros((n_bins, file_cnt))
        for i in xrange(file_cnt):
            # data_fn must be in accordance with the filename given in
            # cluster_odorant_space.py
            data_fn = "affinity_matrix_OR%d_%d.dat" % (n_cluster, i) 
            d = np.loadtxt(folder + "/" +  data_fn)
            d = d.flatten() # we want to count affinities between all OR and all odorants
            global_max = max(global_max, d.max())
            n, bins = np.histogram(d, n_bins)
            all_data[:, i] = n
        
        # 2) average how often a distance has occured in the file_cnt many
        # trials --> create a distribution of pooled distances between ORs and
        # odorants
        for b in xrange(n_bins):
            std_count[b] = all_data[b,:].std()
            mean_count[b] = all_data[b,:].mean()

        # estimating the skewness of the distribution
        alpha_estimate = scipy.stats.skew(mean_count)
        print "Skewness Alpha:", alpha_estimate

        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Affinity distribution for %d ORs" % (n_cluster))
        ax.set_xlabel('Affinities between odorant and ORs')
        ax.set_ylabel('Average number of occurences pooled over %d trials' % file_cnt)
        ax.bar(bins[:-1], mean_count, width=bins[1]-bins[0], yerr=std_count, label="Affinity distribution")

        # fit a gaussian function
        guess_params = [0.2, 0.1, 1.0, .2 * mean_count.max()] # mu, sigma, alpha, amplitude
        bincenters = 0.5*(bins[1:]+bins[:-1])
        x0 = 0
        x1 = n_bins
    #    x1 = int(n_bins / 2.)
        x = bincenters[x0:x1]
        opt_params_skewgauss = leastsq(residuals_skew_normal, guess_params, args=(mean_count[x0:x1], x), maxfev=1000)
    #    opt_params_gauss = leastsq(residuals_gauss, guess_params, args=(mean_count[x0:x1], x), maxfev=1000)
    #    print "opt_mu_sigma", opt_params
        opt_mu = opt_params_skewgauss[0][0]
        opt_sigma = opt_params_skewgauss[0][1]
        opt_alpha = opt_params_skewgauss[0][2]
        opt_max = opt_params_skewgauss[0][3]
        print "Optimal parameters: mu %.2e sigma %.2e alpha %.2e \tmax%.2e" % (opt_mu, opt_sigma, opt_alpha, opt_max)
        opt_skewnormal= opt_max * skew_normal(x, opt_mu, opt_sigma, opt_alpha)
    #    opt_gauss = opt_max * mlab.normpdf(x, opt_mu, opt_sigma)
        ax.plot(x, opt_skewnormal, 'r--', lw=2, label="Skew normal distribution \n(mu=%.1e, sigma=%.1e,\nmax=%.1e, alpha=%.1e" %(opt_mu, opt_sigma, opt_max, opt_alpha))
        pylab.legend()

        ax.set_ylim((0, 1.2 * mean_count.max()))
        fn_out = figure_folder + "affinity_distr_nOR%d_skewfit.png" % (n_cluster)
        print 'Saving figure to:', fn_out
        pylab.savefig(fn_out, dpi=100)



