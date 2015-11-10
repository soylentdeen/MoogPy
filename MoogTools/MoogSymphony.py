import MoogTools
import sys
import AstroUtils
import os
import pyfits
import SpectralTools
import numpy
import random
import string

class Symphony(object ):
    def __init__(self, configFile, moogInstance):
        config = AstroUtils.parse_config(configFile)
        self.flavor = ''.join(random.choice(string.ascii_uppercase) for _ in range(4))
        self.moogInstance = moogInstance

        self.baseconfig = {}
        self.gridconfig = {}
        for key in config.keys():
            if not('grid' in key):
                self.baseconfig[key] = config[key]
            else:
                self.gridconfig[key] = config[key]

        self.generateConfigFiles()

    def generateConfigFiles(self):
        try:
            teff = numpy.array(self.gridconfig['grid_T'].split(','), dtype=numpy.int)
        except:
            teff = [self.baseconfig['Teff']]
        try:
            logg = numpy.array(self.gridconfig['grid_logg'].split(','), dtype=numpy.float)
        except:
            logg = [self.baseconfig['logg']]
        try:
            bfield = numpy.array(self.gridconfig['grid_B'].split(','), dtype=numpy.float)
        except:
            bfield = [self.baseconfig['Bfield']]

        filenames = []
        for temp in teff:
            for grav in logg:
                for b in bfield:
                    configFile = self.baseconfig.copy()
                    configFile["Teff"] = temp
                    configFile["logg"] = grav
                    configFile["Bfield"] = b
                    filename = self.gridconfig["grid_config_directory"]+'Config_T%d_G%.1f_B%.1f_raw.cfg'%(temp, grav, b)
                    AstroUtils.write_config(filename, configFile)
                    filenames.append(filename)
        self.filenames = filenames

    def compose(self):
        path_to_MoogSymphony = os.path.dirname(os.path.abspath(__file__))
        for f in self.filenames:
            os.system('python '+path_to_MoogSymphony+'/generateSpectrum.py '+f+' '+self.flavor+' '+self.moogInstance)
