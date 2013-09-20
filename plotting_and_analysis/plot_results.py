import pylab
import numpy
import sys

if (len(sys.argv) < 2):
    fn = raw_input("Please enter data file to be plotted\n")
else:
    fn = sys.argv[1]

#data = pylab.loadtxt(fn)
data = pylab.loadtxt(fn, skiprows=1)
if (data.ndim == 1):
    x_axis = numpy.arange(data.size)
    pylab.plot(x_axis, data)
else:
    pylab.errorbar(data[:,0], data[:,1], yerr=data[:, 2])
    print 'mean y-value:', data[:, 1].mean()
#    pylab.plot(data[:,0], data[:,1])
#    pylab.scatter(data[:,3], data[:,6])
#    pylab.plot(data[:,3], data[:,6])
pylab.show()
