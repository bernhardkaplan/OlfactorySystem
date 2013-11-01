import os, sys, inspect
# use this if you want to include modules from a subforder
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import numpy as np
import pylab
from FigureCreator import plot_params


noise = np.array([.0, 0.05, .1, .15, .2])
pcompleteness = np.abs(np.array([1., .8, .7, .6, .5, .4 , .3]))
np_noise = np.array([50, 50, 39, 26, 24]) / 50. * 100.
np_complete = np.array([50, 49, 49, 46, 40, 28, 16]) / 50. * 100.

plot_params['figure.subplot.left'] = .17
plot_params['figure.subplot.top'] = .9
pylab.rcParams.update(plot_params)
fig = pylab.figure()
ax = fig.add_subplot(111)

ax.plot(noise, np_noise, '-o', color='b', lw=3)
ax.set_xlabel('Noise in affinity matrix', color='b')
for tl in ax.get_xticklabels():
    tl.set_color('b')

ax2 = ax.twiny()
ax2.plot(pcompleteness, np_complete, '-', marker='^', color='r', lw=3)
ax2.set_xlabel('Pattern completeness', color='r')
for tl in ax2.get_xticklabels():
    tl.set_color('r')

ax.set_ylabel('Percentage of correctly \nrecognized patterns')

output_fn = 'noise_and_completeness_performance.png'
print 'Saving to:', output_fn
pylab.savefig(output_fn, dpi=300)
output_fn = 'noise_and_completeness_performance.pdf'
print 'Saving to:', output_fn
pylab.savefig(output_fn, dpi=300)

pylab.show()