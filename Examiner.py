import Moog960
import MoogTools
import matplotlib.pyplot as pyplot
import SpectralTools

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

vsini = 5.8
R = 40000.0

datafile = './Output/StokesExample_T3600_G4.00_B3.00_raw.fits'
Melody_B3 = Moog960.SyntheticMelody(filename = datafile)
Melody_B3.selectPhrases(selectAll=True)

convolved = Melody_B3.rehearse(vsini=vsini, R=R, returnLabels=True)
spectrum_B3, label = Melody_B3.perform(label=convolved[0])

spectrum_B3.plot(ax=ax)

datafile = './Output/StokesExample_T3600_G4.00_B0.00_raw.fits'
Melody_B0 = Moog960.SyntheticMelody(filename = datafile)
Melody_B0.selectPhrases(selectAll=True)

convolved = Melody_B0.rehearse(vsini=vsini, R=R, returnLabels=True)
spectrum_B0, label = Melody_B0.perform(label=convolved[0])

spectrum_B0.plot(ax=ax)

fig.show()
