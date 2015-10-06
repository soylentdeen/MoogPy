import MoogTools
import matplotlib.pyplot as pyplot


fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#moogPyConfigFile = 'moogPy.cfg'
moogPyConfigFile = 'moogStokesPy.cfg'
Moog = MoogTools.MoogStokes(moogPyConfigFile)
Moog.lineList.writeLineLists(mode="MOOGSTOKES")
Moog.parameterFile.writeParFile()
wavelength, flux = Moog.run(test=False)

ax.plot(wavelength, flux)

fig.show()



