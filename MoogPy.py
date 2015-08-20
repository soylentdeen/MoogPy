import Moog
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import SpectralTools

flux = []
wave = []

def recorder(x, y):
    global flux
    global wave
    wave.append(x)
    flux.append(1.0-y)


configFile = 'gfOptimizer.cfg'
Moog.recorder = recorder

Synth = MoogTools.Configuration(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

Moog.moogsilent()


IM = numpy.zeros((Synth.lineList.numLines, len(wave)))

wavelengths = numpy.array(wave)
nominalSpectrum = numpy.array(flux)

"""

for i in range(Synth.lineList.numLines):
    flux = []
    wave = []
    Synth.lineList.perturbLine(i, 0.3)
    Moog.moogsilent()
    plus = numpy.array(flux)
    flux = []
    wave = []
    Synth.lineList.perturbLine(i, -0.3)
    Moog.moogsilent()
    minus = numpy.array(flux)
    IM[i,:] = plus - minus


U,S,V = scipy.linalg.svd(IM)
"""

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#for i in range(Lines.numLines):
#    ax.plot(wavelengths, IM[i,:])

ax.plot(wavelengths, nominalSpectrum)
ax.plot(Synth.solarSpectrum.wave, Synth.solarSpectrum.flux)

fig.show()

#flux = 1.0 - Moog.linex.d
#wave = Moog.linex.wave1

#ax.plot(wavelengths, 1.0-numpy.array(flux))

#fig.show()

#model.loadIntoFORTRAN(MoogStokes)

"""

Line Parameter Fitting -
generate line list
synthesize Solar
synthesize Arcuturs
for line in line list
    push gf
    synthesize solar/arcturus
    pull gf
    synthesize solar/arcturus
    compute influence function for that line
    append to matrix




Star Paramter Fitting
initial guess

"""
