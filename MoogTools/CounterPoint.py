import astropy.io.fits as pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy

class Player( object ):
    def __init__(self, datafile, ID, resolvingPower):
        self.datafile = datafile
        self.ID = ID
        self.resolvingPower = resolvingPower
        info = pyfits.info(self.datafile, output='')
        self.nSpectra = len(info)-1
        wavestart = []
        wavestop = []
        self.headers = []
        for i in range(self.nSpectra):
            self.headers.append(pyfits.getheader(self.datafile, ext=i+1))
            wavestart.append(self.headers[-1].get('WLSTART'))
            wavestop.append(self.headers[-1].get('WLSTOP'))

        self.wavestart = numpy.array(wavestart)
        self.wavestop = numpy.array(wavestop)
        header = pyfits.getheader(self.datafile)
        self.Teff = header.get("TEFF")
        self.logg = header.get("LOGG")
        self.B = header.get("BFIELD")
        self.vsini = header.get("VSINI")
        self.generate_label()
        self.suppressed = True

    def generate_label(self):
        self.label = "T%4d G%.1f B%.1f V%.1f" % (self.Teff, self.logg, self.B, self.vsini)

    def print_info(self):
        if self.suppressed:
            print("#%d)  %s" % (self.ID, self.label))
        else:
            print("#%d)* %s" % (self.ID, self.label))

    def play(self):
        if self.suppressed:
            return None
        else:
            waves = []
            fluxes = []
            self.loadSpectra()
            # now look for overlaps
            for w, f in zip(self.wave[self.inWlRange==True], self.flux[self.inWlRange==True]):
                waves.append(w)
                fluxes.append(f)
            overlaps = numpy.intersect1d(self.wavestart[self.inWlRange==True], self.wavestop[self.inWlRange==True])
            blue_indices = []
            red_indices = []
            stitched_waves = []
            stitched_fluxes = []
            for overlap in overlaps:
                blue_index = self.wavestop[self.inWlRange==True] == overlap
                red_index = self.wavestart[self.inWlRange==True] == overlap

                for b in numpy.arange(self.nWithinRange)[blue_index]:

                    for r in numpy.arange(self.nWithinRange)[red_index]:
                        stitched_waves.append(numpy.append(waves[b][:-1], waves[r]))
                        stitched_fluxes.append(numpy.append(fluxes[b][:-1], fluxes[r]))
                        if len(blue_indices) == 0:
                            red_indices.append(r)
                    blue_indices.append(b)
            new_waves = []
            new_fluxes = []
            for i in range(self.nWithinRange):
                if not(i in red_indices) and (not(i in blue_indices)):
                    new_waves.append(waves[i])
                    new_fluxes.append(fluxes[i])
            for i in range(len(stitched_waves)):
                new_waves.append(stitched_waves[i])
                new_fluxes.append(stitched_fluxes[i])
            return (self.label, new_waves, new_fluxes)

    def __eq__(self, other):
        try:
            return int(other) == self.ID
        except:
            return self.datafile == other

    def inWlRange(self, wlStart, wlStop):
        inWlRange = []
        for wavestart, wavestop in zip(self.wavestart, self.wavestop):
            if (wavestart < wlStop) & (wavestop > wlStart):
                inWlRange.append(True)
            else:
                inWlRange.append(False)
        self.inWlRange = numpy.array(inWlRange)
        self.nWithinRange = numpy.sum(self.inWlRange)
        if any(self.inWlRange):
            return True
        else:
            return False

    def loadSpectra(self):
        self.wave = []
        self.flux = []
        for i, good in zip(range(self.nSpectra), self.inWlRange):
            if good:
                data = pyfits.getdata(self.datafile, ext=i+1)
                if self.resolvingPower != None:
                    wave, flux = SpectralTools.resample(data.field('Wavelength'), 
                        data.field('Stokes_I'), self.resolvingPower)
                    self.wave.append(wave)
                    self.flux.append(flux)
                else:
                    self.wave.append(data.field('Wavelength'))
                    self.flux.append(data.field('Stokes_I'))
            else:
                self.wave.append([])
                self.flux.append([])
        #self.wave, self.flux = SpectralTools.resample(data[0], data[1], 
        #        self.resolvingPower)
        self.wave = numpy.array(self.wave)
        self.flux = numpy.array(self.flux)

class Orchestra( object ):
    def __init__(self, configFile):
        self.config = AstroUtils.parse_config(configFile)
        self.applyConfigFile()
        self.WhosInThePit()

    def applyConfigFile(self):
        self.watchedDir = self.config["watched_dir"]
        self.wlStart = self.config["wlStart"]
        self.wlStop = self.config["wlStop"]
        if 'observed' in self.config.keys():
            self.resolvingPower = self.config["resolving_power"]
            self.observed = MoogTools.Spectrum(self.config, 'observed')
        else:
            self.observed = None
            self.resolvingPower = None
        self.players = []

    def WhosInThePit(self):
        spectra = glob.glob(self.watchedDir+'*.fits')
        for spectrum in spectra:
            if not(spectrum in self.players):
                player = Player(spectrum, len(self.players), self.resolvingPower)
                if player.inWlRange(self.wlStart, self.wlStop):
                    self.players.append(player)

    def selectPlayers(self):
        self.WhosInThePit()
        for spectrum in self.players:
            spectrum.print_info()

        selection = raw_input("Enter space-separated list of spectra to Toggle :").split()

        for source in selection:
            for spectrum in self.players:
                if spectrum == source:
                    spectrum.suppressed = not(spectrum.suppressed)


    def getEnsemble(self):
        spectra = []
        labels = []
        for player in self.players:
            if not(player.suppressed):
                label, wave, flux = player.play()
                labels.append(label)
                spectra.append([wave, flux])

        return spectra, labels
