import scipy
import numpy

class Atmosphere( object ):
    def __init__(self, df):
        data = open(df, 'r')
        junk = data.readline()
        coords = data.readline().split()
        self.Teff = int(coords[1][2:])
        self.G = float(coords[2][2:])
        self.F = float(coords[3][4:])
        self.m = float(coords[4][3:])
        self.nlayers = int(data.readline().split()[0])
        self.tauref = numpy.zeros(self.nlayers)
        self.T = numpy.zeros(self.nlayers)
        self.Theta = numpy.zeros(self.nlayers)
        self.Tkev = numpy.zeros(self.nlayers)
        self.Tlog = numpy.zeros(self.nlayers)
        self.pgas = numpy.zeros(self.nlayers)
        self.ne = numpy.zeros(self.nlayers)
        self.molweight = numpy.zeros(self.nlayers)
        self.kaprefmass = numpy.zeros(self.nlayers)
        self.rho = numpy.zeros(self.nlayers)
        self.kapref = numpy.zeros(self.nlayers)
        self.mt = 0.0

        for i in range(self.nlayers):
            layer = data.readline().split()
            self.tauref[i] = float(layer[0])
            self.T[i] = float(layer[1])
            self.Theta[i] = 5040./self.T[i]
            self.Tkev[i] = 8.6171e-5*self.T[i]
            self.Tlog[i] = numpy.log10(self.T[i])
            self.pgas[i] = float(layer[2])
            self.ne[i] = float(layer[3])
            self.molweight[i] = float(layer[4])
            self.kaprefmass[i] = float(layer[5])
            self.rho[i] = self.pgas[i]*self.molweight[i]*1.6606e-24/(1.38054e-16*self.T[i])
            self.kapref[i] = self.kaprefmass[i] * self.rho[i]
        self.mt = float(data.readline().split()[0])
        data.close()

    def loadIntoFORTRAN(self, FortObj):
        for i in range(self.nlayers):
            FortObj.atmos.tauref[i] = self.tauref[i]
            FortObj.atmos.t[i] = self.T[i]
            FortObj.atmos.theta[i] = self.Theta[i]
            FortObj.atmos.tkev[i] = self.Tkev[i]
            FortObj.atmos.tlog[i] = self.Tlog[i]
            FortObj.atmos.pgas[i] = self.pgas[i]
            FortObj.atmos.ne[i] = self.ne[i]
            FortObj.atmos.molweight[i] = self.molweight[i]
            FortObj.atmos.kapref[i] = self.kapref[i]


class Params( object ):
    def __init__(self, parfile):
        self.filename = parfile

class AtomicLines( object ):
    def __init__(self, df):
        self.filename = df

class MolecularLines( object ):
    def __init__(self, df):
        self.filename = df
