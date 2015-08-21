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
Synth.solarSpectrum.fixDiscontinuities()
Synth.solarSpectrum.flipWavelength()

solarWave = Synth.solarSpectrum.wave

#"""
Moog.moogsilent()

#
# Construction of IM:
# log-gfs for each line - check
# Continuum Level - check
# WL shift - check
# Instrumental Broadening
IM = numpy.zeros((Synth.lineList.numLines+2, len(solarWave)))

wavelengths = numpy.array(wave)
nominalSpectrum = numpy.array(flux)

newNominalSpectrum = SpectralTools.binSpectrum(nominalSpectrum, wavelengths,
        solarWave)
#"""
#"""

for i in range(Synth.lineList.numLines):
    flux = []
    wave = []
    Synth.lineList.perturbLine(i, 0.3)
    Moog.moogsilent()
    plus = SpectralTools.binSpectrum(numpy.array(flux), wavelengths, solarWave)
    flux = []
    wave = []
    Synth.lineList.perturbLine(i, -0.3)
    Moog.moogsilent()
    minus = SpectralTools.binSpectrum(numpy.array(flux), wavelengths, solarWave)
    IM[i,:] = plus - minus

#Continuum Level
plus = newNominalSpectrum.copy()+ 0.005
mins = newNominalSpectrum.copy()- 0.005
IM[Synth.lineList.numLines, :] = plus-minus

plus = SpectralTools.binSpectrum(nominalSpectrum, wavelengths+0.1, solarWave)
minus = SpectralTools.binSpectrum(nominalSpectrum, wavelengths-0.1, solarWave)
edges = (plus !=0) & (minus != 0)
IM[Synth.lineList.numLines+1, edges] = plus[edges]-minus[edges]

#U,S,V = scipy.linalg.svd(IM)
#"""

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

"""
ax.plot(wavelengths, nominalSpectrum)
ax.plot(Synth.solarSpectrum.wave, newNominalSpectrum)
ax.plot(Synth.solarSpectrum.wave, Synth.solarSpectrum.flux)
#"""

#"""
for i in range(Synth.lineList.numLines+2):
    ax.plot(solarWave, IM[i,:])
#"""
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
