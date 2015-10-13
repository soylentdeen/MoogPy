import MoogTools
import matplotlib.pyplot as pyplot
import numpy
import pyfits


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

moogPyConfigFile = 'TraceStokes.cfg'
Moog = MoogTools.MoogStokesSpectrum(moogPyConfigFile, 'trace')
tau, I, Q, U, V, cont = Moog.trace(save=True)

ax.plot(tau, numpy.array(I)/numpy.array(cont))

fig.show()



