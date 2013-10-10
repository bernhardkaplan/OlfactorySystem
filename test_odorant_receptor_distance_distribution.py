import numpy as np
import matplotlib.mlab as mlab
import pylab

def distance_generator():

#    p = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
#         1.909487e+01, 1.110809e+01, 3.353855e+00, \
#         4.188897e+00, 4.088460e+01, 4.966478e-01]

    p = [162.310869565, 6.67080434783, 1.98630434783,\
          19.8056521739, 10.8089130435, 3.32682608696, \
           4.4382173913, 40.8932608696, 0.456293478261] # these values are taken from clustering the odorant space with 40 ORs
    w1 = p[0]
    mu1 = p[1]
    sigma1 = p[2]
    w2 = p[3]
    mu2 = p[4]
    sigma2 = p[5]
    w3 = p[6]
    mu3 = p[7]
    sigma3 = p[8]
    which_gauss = np.random.uniform(0, 1.)
    if which_gauss < p1:
#        print 'G1 ', 
        return max(0, np.random.normal(mu1, sigma1))
    elif (which_gauss < p2 + p1):
#        print 'G2 ', 
        return max(0, np.random.normal(mu2, sigma2))
    elif (which_gauss < p3 + p2 + p1):
#        print 'G3 ', 
        return max(0, np.random.normal(mu3, sigma3))
    else:
        print '\nERROR\n'


def get_prob_for_distributions():
    """
    Based on the integral of the three normal distributions,
    the likelihood from which of the three distributions a distance is to be drawn 
    is calculated here.
    Returns the three probabilities for the three distributions.
    """
    p = [1.631787e+02, 6.670855e+00, 1.977871e+00, \
         1.909487e+01, 1.110809e+01, 3.353855e+00, \
         4.188897e+00, 4.088460e+01, 4.966478e-01]
    w1 = p[0]
    mu1 = p[1]
    sigma1 = p[2]
    w2 = p[3]
    mu2 = p[4]
    sigma2 = p[5]
    w3 = p[6]
    mu3 = p[7]
    sigma3 = p[8]

    dist_range = (0, 4.330310991999920844e+01)
    x = np.linspace(dist_range[0], dist_range[1], 1000)
    A1 = np.array(w1 * mlab.normpdf(x, mu1, sigma1)).sum()
    A2 = np.array(w2 * mlab.normpdf(x, mu2, sigma2)).sum()
    A3 = np.array(w3 * mlab.normpdf(x, mu3, sigma3)).sum()
    print 'debug A', A1, A2, A3
    p1 = A1 / (A1 + A2 + A3)
    p2 = A2 / (A1 + A2 + A3)
    p3 = A3 / (A1 + A2 + A3)
    return p1, p2, p3

p1, p2, p3 = get_prob_for_distributions()
print 'Debug p1 p2 p3', p1, p2, p3, p1 + p2 + p3
np.random.seed(0)

n_or = 500 
n_pattern = 1000
n_samples = n_or * n_pattern
samples = np.zeros(n_samples)

for i_ in xrange(n_samples):
    samples[i_] = distance_generator()
count, bins = np.histogram(samples, bins=200)

count_norm = count / float(n_samples)
expected_value = np.sum(count_norm * bins[:-1])
pylab.bar(bins[:-1], count_norm, width=bins[1]-bins[0])
print 'debug counts', count, n_samples
print 'debug normalized counts', count_norm
print 'bins * count_norm', bins[:-1] * count_norm, (bins[:-1] * count_norm).sum()
print 'Expected value:', expected_value
#pylab.bar(bins[:-1], count, width=bins[1]-bins[0])


pylab.show()

