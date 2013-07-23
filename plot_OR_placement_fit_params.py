
import pylab
import numpy as np
import sys

if (len(sys.argv) < 2):
    fn = raw_input("Please enter data file to be plotted\n")
else:
    fn = sys.argv[1]

rcp = {'figure.subplot.wspace': 0.40, 'figure.subplot.hspace':.4}
pylab.rcParams.update(rcp)
fig = pylab.figure()
ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)
data = pylab.loadtxt(fn)

ax1.scatter(data[:, 0], data[:, 1])
ax2.scatter(data[:, 0], data[:, 2])
ax3.scatter(data[:, 0], data[:, 3])
ax4.scatter(data[:, 0], data[:, 4])
ax5.scatter(data[:, 0], data[:, 5])
ax6.scatter(data[:, 0], data[:, 6])
ax7.scatter(data[:, 0], data[:, 7])
ax8.scatter(data[:, 0], data[:, 8])
ax9.scatter(data[:, 0], data[:, 9])




pylab.show()


