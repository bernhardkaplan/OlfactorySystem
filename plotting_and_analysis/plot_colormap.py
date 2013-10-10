import pylab
import numpy as np
import sys
import os

def get_mean_max_min_median(d):
    print 'In row   mean    max     min     mean of smallest 10 percent   median'
    for i in xrange(d[:, 0].size):
        idx_sorted = d[i, :].argsort()
        lowest = d[i, idx_sorted[:10]]
        print i, d[i, :].mean(), d[i, :].max(), d[i, :].min(), lowest.mean(), np.median(d[i, :])
        

def get_units_above_thresh(d, thresh):

    n_above_thresh = np.zeros(d[:, 0].size)
    for i in xrange(d[:, 0].size):
        above_thresh = (d[i, :] > thresh).nonzero()[0]
        n_above_thresh[i] = above_thresh.size
        print 'In row %d: %d units have value higher than %.2e' % (i, n_above_thresh[i], thresh), above_thresh
    print 'On average, per row %.2f +- %.2f units have higher values than %.2e, i.e. %.1f +- %.1f percent' % (n_above_thresh.mean(), n_above_thresh.std(), thresh, \
            n_above_thresh.mean() / d[0, :].size * 100., n_above_thresh.std() / d[0, :].size * 100.)
    print 'In total %d units have higher values than %.2e' % (n_above_thresh.sum(), thresh)


if (len(sys.argv) < 2):
    print("Please enter data folder to be analysed after the script name")
    print("e.g.\npython analyse_data_folder.py data_today/ ")
    exit(1)
else:
    fn = sys.argv[1]



path = os.path.abspath(fn)

print "loading data ...."

# if it's only one line:
#d = np.loadtxt(path)
#data = np.zeros((1, d.size))
#for i in xrange(d.size):
#    data[0, i] = d[i]

try:
    data = np.loadtxt(path, delimiter=",")#.transpose()
except:
    data = np.loadtxt(path)

get_mean_max_min_median(data)
thresh = 0.0001
get_units_above_thresh(data, thresh)

# if you want to take the log of the data to plot
#n_row = data[:, 0].size
#n_col = data[0, :].size
#log_data = np.zeros((n_row, n_col))
#for i in xrange(n_row):
#    for j in xrange(n_col):
#        if data[i, j] > 0:
#            log_data[i, j] = np.log(data[i, j])

#data_rev = np.zeros(data.shape)
#n_row = data[:, 0].size - 1
#for row in xrange(data[:, 0].size):
#    data_rev[n_row - row, :] = data[row, :]

fig = pylab.figure()
ax = fig.add_subplot(111)
print "plotting ...."
#cax = ax.imshow(data[:,:12])
#cax = ax.pcolor(data, edgecolor='k', linewidths='1')

cax = ax.pcolormesh(data)#, edgecolor='k', linewidths='1')
#cax = ax.pcolor(data, cmap='binary')
#cax = ax.pcolor(data, cmap='RdBu')

#cax = ax.pcolor(log_data)#, edgecolor='k', linewidths='1')


pylab.ylim(0, data.shape[0])
pylab.xlim(0, data.shape[1])
pylab.colorbar(cax)

#plot_fn = "testfig.png"
#print "saving ....", plot_fn
#pylab.savefig(plot_fn)

pylab.show()
