import pylab
import numpy
import sys
import os

num_files = len(sys.argv) - 1

if (num_files == 0):
    fn = raw_input("Please enter data file to be plotted\n")

#print sys.argv[1:]

pylab.figure()
for fn in sys.argv[1:]:
    data = pylab.loadtxt(fn, skiprows=1)


    if (data.ndim == 1):
        x_axis = numpy.arange(data.size)
        pylab.plot(x_axis, data, label=fn)
#        pylab.plot(x_axis, data)
    else:
        pylab.plot(data[:,0], data[:,1], lw=2, label=fn)#, yerr=data[:,3])
        #pylab.title("%s for n_rows=1000, averaged over all processes" % (title))
        pylab.xticks(size='large')
        pylab.yticks(size='large')
        #pylab.savefig("%s_nrows1000.eps" % (title))
    print fn

#pylab.legend()
pylab.show()
