import MoogTools
import matplotlib.pyplot as pyplot


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

moogPyConfigFile = 'StokesExample.cfg'
Moog = MoogTools.MoogStokes(moogPyConfigFile, fileBase = 'example')
Moog.run(saveRaw=True)


ax.plot(wavelength, flux)

fig.show()



