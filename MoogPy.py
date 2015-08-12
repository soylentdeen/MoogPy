import Moog
import PyMoog
import numpy
import scipy
import matplotlib.pyplot as pyplot

flux = []
wavelengths = []

def recorder(x, y):
    global flux
    global wavelengths
    wavelengths.append(x)
    flux.append(y)


#datadir = '/home/deen/Data/MoogStokes/Atmospheres/MARCS'
#modelfile = datadir+'/MARCS_T3400_G4.0_M0.0_t2.0.md'

#paramFile = 'ParFile'
#atomicFile = 'AtomicFile'
#molecFile = 'MolecFile'

#model = PyMoog.Atmosphere(modelfile)

#params = PyMoog.Params(paramFile)

#atomicLines = PyMoog.AtomicLines(atomicFile)
#moleucularLines = PyMoog.MolecularLines(molecFile)


Moog.recorder = recorder

Moog.moogsilent()

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#flux = 1.0 - Moog.linex.d
#wave = Moog.linex.wave1

ax.plot(wavelengths, 1.0-numpy.array(flux))

fig.show()

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
