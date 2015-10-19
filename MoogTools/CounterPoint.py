import pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools

class Player( object ):
    def __init__(self, datafile, ID, resolvingPower):
        self.datafile = datafile
        self.ID = ID
        self.resolvingPower = resolvingPower
        data = pyfits.getdata(self.datafile)
        self.wave, self.flux = SpectralTools.resample(data[0], data[1], 
                self.resolvingPower)
        header = pyfits.getheader(self.datafile)
        self.Teff = header.get("TEFF")
        self.logg = header.get("LOGG")
        self.B = header.get("BFIELD")
        self.vsini = header.get("VSINI")
        self.generate_label()
        self.suppressed = False

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
            return (self.label, self.wave, self.flux)

    def __eq__(self, other):
        try:
            return int(other) == self.ID
        except:
            return self.datafile == other

    def inWlRange(self, wlStart, wlStop):
        return (self.wave[0] < wlStop) & (self.wave[-1] > wlStart)

class Orchestra( object ):
    def __init__(self, configFile):
        self.config = AstroUtils.parse_config(configFile)
        self.applyConfigFile()
        self.WhosInThePit()

    def applyConfigFile(self):
        self.watchedDir = self.config["watched_dir"]
        self.wlStart = self.config["wlStart"]
        self.wlStop = self.config["wlStop"]
        self.resolvingPower = self.config["resolving_power"]
        self.observed = MoogTools.Spectrum(self.config, 'observed')
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

        selection = raw_input("Enter space-separated list of spectra to Toggle :")

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
