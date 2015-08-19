import Moog
import PyMoog
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools

flux = []
wave = []

def recorder(x, y):
    global flux
    global wave
    wave.append(x)
    flux.append(y)
    print wave[-1], flux[-1]
    if flux[-1] < 0:
        raw_input()


configFile = 'gfOptimizer.cfg'
Moog.recorder = recorder
Lines = MoogTools.LineList(configFile, 11700, 11730, 0.0, molecules=False)
Lines.writeLineLists()
Moog.moogsilent()

IM = numpy.zeros((Lines.numLines, len(wave)))

wavelengths = numpy.array(wave)
nominalSpectrum = numpy.array(flux)

for i in range(Lines.numLines):
    flux = []
    wave = []
    #Moog.rewind()
    Lines.perturbLine(i, 0.3)
    Moog.moogsilent()
    plus = numpy.array(flux)
    flux = []
    wave = []
    #Moog.rewind()
    Lines.perturbLine(i, -0.3)
    Moog.moogsilent()
    minus = numpy.array(flux)
    IM[i,:] = plus - minus


U,S,V = scipy.linalg.svd(IM)

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#for i in range(Lines.numLines):
#    ax.plot(wavelengths, IM[i,:])

ax.plot(wavelengths, nominalSpectrum)

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
