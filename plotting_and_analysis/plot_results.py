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
    pylab.scatter(data[:,3], data[:,6])
#    pylab.plot(data[:,3], data[:,6])
pylab.show()
