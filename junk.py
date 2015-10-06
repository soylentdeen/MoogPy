import scipy
import numpy
import matplotlib.pyplot as pyplot

df = '/home/deen/Code/Python/MoogPy/data/AbsorptionLines/VALD/vald_master.txt'

data = open(df, 'r')

wl = []
sf = []
hbar = 4.135667662e-15
c = 2.99792458e18

for line in data.readlines():
    try:
        parsed = line.split(',')
        wave = numpy.float(parsed[1])/10000.0
        calc_wl = c/((numpy.float(parsed[3]) - numpy.float(parsed[5]))/hbar)
        print asdf
        sf.append(calc_wl/wave)
        wl.append(wave)
        #if (1.9 < wave) & (2.5 > wave):
        #    wl.append(wave)
    except:
        pass

wl = numpy.array(wl)

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

ax.plot(wl, sf)
#ax.hist(wl, bins = 100)

fig.show()

