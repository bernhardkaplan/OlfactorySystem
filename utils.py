import numpy
import numpy.random as rnd
import math
import matplotlib.mlab as mlab
from scipy.special import erf
from scipy.spatial import distance

"""
functions_math:
    This module defines some helper functions.
    - normal distribution (x, mu, sigma)
"""
def euclidean(x, y):
    return distance.euclidean(x, y)
    
def normal_distr(x, mu, sigma, A=1.):
    val = A / numpy.sqrt(2 * numpy.pi * float(sigma)**2) * numpy.exp((-1) * (x - mu)**2 / (2 * float(sigma)**2) )
    return val

def pdf(x):
    return 1/numpy.sqrt(2*numpy.pi) * numpy.exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/numpy.sqrt(2))) / 2

def skew_distr(mu, sigma, alpha):
    #for reference, see: http://azzalini.stat.unipd.it/SN/faq-r.html
    # DOI:  10.1111/1467-9868.00194
    # or Statistical applications of the multivariate skew-normal distribution, Adelchi Azzalini, Antonella Capitanio, J. Royal Statistical Society, series B, 61 (1999) 579-602
    u0 = rnd.normal(0, 1)
    v = rnd.normal(0, 1)
    delta = alpha / numpy.sqrt(1 + alpha**2)
    u1 = delta * u0 + numpy.sqrt(1 - delta**2) * v
    if ( u0 >= 0.):
        return u1 * sigma + mu
    else:
        return -u1 * sigma + mu

def skew_normal(x, mu, sigma, alpha):
    # alpha  = skewness parameter
    return (1. / sigma) * mlab.normpdf((x - mu)/sigma, 0., 1.) * (1. + erf(alpha * (x - mu) / (sigma * numpy.sqrt(2))) )

def get_trimodal_normal(w1, mu1, sigma1, w2, mu2, sigma2, w3, mu3, sigma3):
    x = rnd.uniform()
    if (x < w1):
        return rnd.normal(mu1, sigma1)
    elif (x < w1 + w2):
        return rnd.normal(mu2, sigma2)
    else:
        return rnd.normal(mu3, sigma3)


    return w1 * mlab.normpdf(x, mu1, sigma1) + w2 * mlab.normpdf(x, mu2, sigma2) + w3 * mlab.normpdf(x, mu3, sigma3)

def gor_func(x, a1, a2, exp):
    """this function calculates the gor for  a given x """
    return a1 * x**exp + a2

def draw_connection(p, w, noise=0):
    """
    Decide whether a connection is drawn, given the possibility p.
    w : is the mean weight for the type of connection to be drawn.
    noise : is an optional argument and stands for the relative sigma of the normal distributino
    i.e. noise = 0.1 means sigma = 0.1 * w
    """
    if (p > rnd.random()):
        if (noise != 0 and noise > 0):
            weight = rnd.normal(w, noise)
            # check if sign of weight changed
            # if so, return 0
            if (numpy.sign(weight) != numpy.sign(w)):
                return 0
            return weight
        elif (noise < 0): # stupid user, noise should be > 0
            print "WARNING, negative noise given to functions.draw_connection(p, w, noise)!"
            noise *= (-1.0)
            weight = rnd.normal(w, noise)
            if (numpy.sign(weight) != numpy.sign(w)):
                return 0
            return weight
        elif (noise == 0):
            return w
    return 0

def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (math.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

def get_p(l):
    l = list(l)
    p_list = []
    for i in l:
        p = float(l.count(i)) / len(l)
        p_list.append(p)
    return p_list

def shannon_entropy(l):
    """Return the Shannon entropy of random variable with probability
    vector l."""
    return sum([-p*numpy.log2(p) for p in l if p > 0])

def get_Q_wta(l, alpha=2):
    """
    Return the WTA characteristics of an output rate vector.
    The return value is equivalent to Q_ci defined in B.Kaplan's diploma thesis: p. 42 ff
    The smaller alpha is, the more strict the assessment of the spike output is with respect to a perfect winner-take-
    all output. The larger alpha is, the smaller s* has to be for a Qci > 0: e.g. alpha = 5 requires
    that the most active neuron fires more that 20 % of the total spike output in order that
    Qci > 0.
    """
    assert (alpha > 1), "get_Q_wta: alpha must be > 1!"
    l = numpy.array(l)
    s_max = l.max()
    s_sum = l.sum()
    print "debug", s_max, s_sum, s_max / s_sum
    if (s_max == s_sum):
        return 0
    return ((alpha * s_max - s_sum) / s_sum) * (1 / (float(alpha) - 1)) * (s_sum / s_max)


def gauss(x, mu, sigma):
    return numpy.exp( - (x - mu)**2 / (2 * sigma ** 2))


# ------ for plotting -------------
def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (numpy.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size


def get_figsize_A4():
    fig_width = 8.27
    fig_height = 11.69
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

def get_figsize_A4_landscape():
    fig_width = 11.69 
    fig_height = 8.27
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size


def get_figsize_A5_landscape():
    fig_width = 8.27
    fig_height = 5.83
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size


def distribute_numbers(N, n_proc, pid):
    """
    distribute N elements (cells, whatever...) among n_proc processes
    returns tuple (limit_min, limit_max) for the process with pid
    """
    # ------ divide the N elements among the processes: each process deals with a subset
    n_mine = int(N / n_proc)
    R = N % n_proc
    offset = min(pid, R)   # possibly some processes have to deal with more elements
    limit_min = int(pid * n_mine + offset)
    if (pid < R):
        limit_max = int(limit_min + n_mine + 1)
    else:
        limit_max = int(limit_min + n_mine)
    return (limit_min, limit_max)


def distribute_list(l, n_proc, pid):
    """
    This function picks a number of elements from the list l so that all of the n_proc processes get a different list.
    Arguments:
        l = list of numbers to be distributed
        n_proc = number of processors
        pid = the target processor 
    """
    N = len(l)
    (limit_min, limit_max) = distribute_numbers(N, n_proc, pid) # this gives the indices which elements to pick
    my_list = []
    for i in xrange(limit_min, limit_max):
        my_list.append(l[i])
    return my_list


def orn_sigmoid(x, a, b, c, d, tau):
    """
    function which generates a sigmoidal dose-response curve
    """
    # d = limit value for x -> - infinity
    # a, b = limit value for x -> + infinity
    # tau, c = position for transition
    f_x = a / (b + d * numpy.exp(-tau * (x - c)))
    return f_x


def get_sources(filename, cell_id):
    d = numpy.loadtxt(filename, skiprows=1)
    srcs = []
    for row in xrange(d[:,0].size):
        tgt = d[row, 1]
        if (tgt == cell_id):
            srcs.append(d[row, 0])
    return srcs


def get_targets(filename, cell_id):
    d = numpy.loadtxt(filename, skiprows=1)
    tgts = []
    for row in xrange(d[:,0].size):
        src = d[row, 0]
        if (src == cell_id):
            tgts.append(d[row, 1])
    return tgts 


def fill_up_nspike_data(fn, n_cells, offset):
    data = numpy.loadtxt(fn) 
    data_new = numpy.zeros((n_cells, data.shape[1]))
    for i in xrange(data.shape[0]):
        gid = data[i, 0]
        index = gid - offset
        data_new[index, 1] = data[i, 1]
    return data_new

def wta(d):
    """
    d : 1-D numpy.array e.g. of activity
    returns the array after the winner-take-all operation
    """
    winner = numpy.argmax(d)
    d_new = numpy.zeros(d.size)
    d_new[winner] = 1.0
    return d_new

    

def get_rnd_targets(gid_max, n):
    """
    gid_max: upper boundary
    n: numbe of random integers to be drawn within range [0, gid_max]
    """
    tgts = rnd.randint(0, gid_max, n)
    tgts = numpy.unique(tgts).tolist()
    while len(tgts) < n:
        # check if rnd_int is already in l 
        tgts.append(rnd.randint(0, gid_max))
        tgts = numpy.unique(tgts).tolist()
    return tgts

