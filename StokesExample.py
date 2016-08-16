import MoogTools
import Moog960
import SpectralTools
import astropy.io.fits as pyfits
import matplotlib.pyplot as pyplot


fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
#ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])

moogPyConfigFile = 'StokesExample.cfg'
MoogRaw = MoogTools.MoogStokes(moogPyConfigFile, fileBase = 'example')
MoogRaw.run(saveRaw=True)



MoogRaw.Phrase.rawData[0].plot(ax=ax1)
MoogRaw.Phrase.rawData[1].plot(ax=ax1)
MoogRaw.Phrase.rawData[2].plot(ax=ax1)
MoogRaw.Phrase.rawData[3].plot(ax=ax1)
MoogRaw.Phrase.rawData[4].plot(ax=ax1)
MoogRaw.Phrase.rawData[5].plot(ax=ax1)
MoogRaw.Phrase.rawData[6].plot(ax=ax1)
#Moog.Phrase.convolvedData[0].plot(ax=ax2)

#MoogProcessed = Moog960.SyntheticPhrase.fromFile(None, filename='./Output/StokesExample_T3600_G4.00_B0.00_raw.fits')

fig.show()



