import pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy
import os
import matplotlib.lines as Lines

class Phrase( object ):
    def __init__(self, rawData=None, diskInt = 'BEACHBALL'):
        self.rawData = rawData
        self.wlStart = rawData[0].header.get('WLSTART')
        self.wlStop = rawData[0].header.get('WLSTOP')
        if diskInt == 'BEACHBALL':
            self.processedData = SpectralTools.BeachBall(parent=self)
        elif diskInt == 'DISKOBALL':
            self.processedData = SpectralTools.DiskoBall(parent=self)

    @classmethod
    def fromFile(self, hdr, data=None, filename=None, ext=None):
        rawData = []
        rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
            filename=filename, ext=ext))
        return self(rawData=rawData)

    def addRawSpectrum(self, hdr, data=None, filename=None, ext=None):
        self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
            filename=filename, ext=ext))

    def owns(self, hdr):
        if ((self.wlStart == hdr.get('WLSTART')) & 
                (self.wlStop == hdr.get("WLSTOP"))):
            return True
        return False

    def inWlRange(self, wlStart, wlStop):
        return ((self.wlStart < wlStop) & (self.wlStop > wlStart))

    #def integrate(self, vsini=0.0):
    #    self.processedData.diskInt(vsini=vsini)

    def rehearse(self, vsini=0.0, R=0):
        self.processedData.resample(vsini=vsini, R=R)

    def perform(self, vsini= 0.0, R = 0.0):
        return self.processedData.yank(vsini=vsini, R = R)

    def saveRaw(self, filename = None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.rawData:
            hdr = spectrum.header.copy()
            SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
            SpectrumHDU.name = "%.4fA - %.4fA PHI=%.3f MU=%.3f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get("PHI_ANGLE"), hdr.get("MU"))

            HDUs.append(SpectrumHDU)

        if filename == None:
            return HDUs

        if os.path.exists(filename):    #file exists!  get ready to append!
            while os.path.exists(filename+'.lock'):
                print("Gnarly dude!  The file is locked!  I'll just hang out here for a while and wait")
                time.sleep(0.1)
            with open(filename+'.lock', 'w'):
                os.utime(filename+'.lock', None)
            HDUList = pyfits.open(filename, mode='update')
            for spectrum in HDUs:
                try:
                    HDUList.pop(HDUList.index_of(spectrum.name))
                except:
                    pass
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.close()
            os.remove(filename+'.lock')
        else:
            HDUList = pyfits.HDUList()
            primary = pyfits.PrimaryHDU()
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)


    def saveInterpolated(self, filename = None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.processedData.interpolated:
            hdr = spectrum.header.copy()
            SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
            SpectrumHDU.name = "%.4fA - %.4fA PHI=%.3f MU=%.3f DELTAV=%.3f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get("PHI_ANGLE"), hdr.get("MU"), hdr.get('DELTAV'))

            HDUs.append(SpectrumHDU)

        if filename == None:
            return HDUs

        if os.path.exists(filename):    #file exists!  get ready to append!
            while os.path.exists(filename+'.lock'):
                print("Gnarly dude!  The file is locked!  I'll just hang out here for a while and wait")
                time.sleep(0.1)
            with open(filename+'.lock', 'w'):
                os.utime(filename+'.lock', None)
            HDUList = pyfits.open(filename, mode='update')
            for spectrum in HDUs:
                try:
                    HDUList.pop(HDUList.index_of(spectrum.name))
                except:
                    pass
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.close()
            os.remove(filename+'.lock')
        else:
            HDUList = pyfits.HDUList()
            primary = pyfits.PrimaryHDU()
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)

    def saveIntegrated(self, filename = None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.processedData.integrated:
            hdr = spectrum.header.copy()
            SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
            SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'))

            HDUs.append(SpectrumHDU)

        if filename == None:
            return HDUs

        if os.path.exists(filename):    #file exists!  get ready to append!
            while os.path.exists(filename+'.lock'):
                print("Gnarly dude!  The file is locked!  I'll just hang out here for a while and wait")
                time.sleep(0.1)
            with open(filename+'.lock', 'w'):
                os.utime(filename+'.lock', None)
            HDUList = pyfits.open(filename, mode='update')
            for spectrum in HDUs:
                try:
                    HDUList.pop(HDUList.index_of(spectrum.name))
                except:
                    pass
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.close()
            os.remove(filename+'.lock')
        else:
            HDUList = pyfits.HDUList()
            primary = pyfits.PrimaryHDU()
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)

    def saveConvolved(self, filename = None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.processedData.convolved:
            hdr = spectrum.header.copy()
            SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
            SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f R=%.1f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'), hdr.get('RESOLVING_POWER'))

            HDUs.append(SpectrumHDU)

        if filename == None:
            return HDUs

        if os.path.exists(filename):    #file exists!  get ready to append!
            while os.path.exists(filename+'.lock'):
                print("Gnarly dude!  The file is locked!  I'll just hang out here for a while and wait")
                time.sleep(0.1)
            with open(filename+'.lock', 'w'):
                os.utime(filename+'.lock', None)
            HDUList = pyfits.open(filename, mode='update')
            for spectrum in HDUs:
                try:
                    HDUList.pop(HDUList.index_of(spectrum.name))
                except:
                    pass
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.close()
            os.remove(filename+'.lock')
        else:
            HDUList = pyfits.HDUList()
            primary = pyfits.PrimaryHDU()
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)

class Melody( object ):
    def __init__(self, phrases = [], filename=None):
        self.phrases = []
        self.filename = filename
        self.loadMelody()
        self.muted = True
        self.nPhrases = len(self.phrases)

    def loadMelody(self):
        info = pyfits.info(self.filename, output='')
        self.nSpectra = len(info)-1
        self.header = pyfits.getheader(self.filename, ext=0)
        self.Teff = self.header.get("TEFF")
        self.logg = self.header.get("LOGG")
        self.B = self.header.get("BFIELD")

        for i in range(self.nSpectra):
            added = False
            hdr = pyfits.getheader(self.filename, ext=i+1)
            #data = pyfits.getdata(self.filename, ext=i+1)
            for phrase in self.phrases:
                if phrase.owns(hdr):
                    phrase.addRawSpectrum(hdr, data=None, filename=self.filename,
                            ext=i+1)
                    #phrase.addRawSpectrum(hdr, data)
                    added=True
                    break
            if not(added):
                self.phrases.append(Phrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1))

        #for phrase in self.phrases:
        #    phrase.processedData.loadData()

    def addPhrase(self, phrases):
        for phrase in phrases:
            self.phrases.append(phrase)
        self.nPhrases = len(self.phrases)

    def selectPhrases(self, wlRange=[]):
        self.selectedPhrases = []
        for phrase in self.phrases:
            self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                wlStop=wlRange[1]))
        print("T=%d G=%d")

    def inParameterRange(self, TeffRange=[], loggRange=[], BfieldRange=[]):
        self.muted = False
        try:
            if ((self.Teff < TeffRange[0]) or (self.Teff > TeffRange[1])):
                self.muted = True
        except:
            pass
        try:
            if ((self.logg < loggRange[0]) | (self.logg > loggRange[1])):
                self.muted = True
        except:
            pass
        try:
            if ((self.B < BfieldRange[0]) | (self.B > BfieldRange[1])):
                self.muted = True
        except:
            pass

    def rehearse(self, vsini = 0.0, R = 0):
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].rehearse(vsini = vsini, R=R)

    def perform(self, vsini = 0.0, R = 0.0):
        spectra = []
        label = "Teff = %d K log g = %.2f Bfield = %.2f kG" % (self.Teff, self.logg, self.B)
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                spectra.append(self.phrases[i].perform(vsini=vsini, R=R))


        return spectra, label

    def record(self, filename=''):
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].record(filename)

class Score( object ):
    """
    
    """
    def __init__(self, melodies = [], directory=None):
        self.melodies = melodies
        self.directory = directory
        self.loadMelodies()

    def loadMelodies(self):
        melodyFiles = glob.glob(self.directory+'*raw.fits')
        for melody in melodyFiles:
            print("%s" % melody)
            self.melodies.append(Melody(filename=melody))

    def setWlRange(self, wlStart, wlStop):
        for melody in self.melodies:
            melody.selectPhrases(wlStart, wlStop)

    def getMelodyParams(self):
        Teff = []
        logg = []
        B = []
        for melody in self.melodies:
            Teff.append(melody.Teff)
            logg.append(melody.logg)
            B.append(melody.B)

        return numpy.unique(Teff), numpy.unique(logg), numpy.unique(B)


    def selectEnsemble(self, T=[], G=[], B=[]):
        print("%s, %s %s" % (T, G, B))
        for melody in self.melodies:
            if (melody.Teff in T) and (melody.logg in G) and (melody.B in B):
                melody.muted=False
                print("Loud: %d, %.1f, %.1f" % (melody.Teff, melody.logg, melody.B))
            else:
                melody.muted=True
                print("Mute: %d, %.1f, %.1f" % (melody.Teff, melody.logg, melody.B))

    def selectMelodies(self, TeffRange = [], loggRange = [], BfieldRange=[],
            wlRange=[]):
        for melody in self.melodies:
            melody.inParameterRange(TeffRange=TeffRange,
                loggRange=loggRange, BfieldRange=BfieldRange)
            if not(melody.muted):
                melody.selectPhrases(wlRange=wlRange)

    def rehearse(self, vsini=0.0, R=0.0):
        '''
        Score.rehearse(vsini=0.0, R=0.0) - For the melodies and phrases selected by
             the parameter and wavelength ranges, generate spectra corresponding to
             the requested vsini and resolving power.

             INPUT:
                 vsini - rotational broadening - km/s
                 R - resolving power

             OUTPUT:
                 none

        The 
        '''
        for melody in self.melodies:
            if not(melody.muted):
                melody.rehearse(vsini=vsini, R=R)

    def perform(self, vsini = 0.0, R = 0.0):
        spectra = []
        labels = []
        for melody in self.melodies:
            if not(melody.muted):
                sp, label = melody.perform(vsini=vsini, R=R)
                spectra.append(sp)
                labels.append(label)

        return spectra, labels

    def record(self, vsini = 0.0, R = 0.0, filename=''):
        for melody in self.melodies:
            if not(melody.muted):
                melody.record(outfilename)

class Moog960( object ):
    def __init__(self, configFile):
        self.config = AstroUtils.parse_config(configFile)
        self.applyConfigFile()

    def applyConfigFile(self):
        self.watchedDir = self.config["watched_dir"]
        self.wlRange = numpy.array(self.config["wlRange"].split(), dtype=numpy.int)
        keys = self.config.keys()
        if "TeffRange" in keys:
            self.TeffRange = numpy.array(self.config["TeffRange"].split(),
                    dtype=numpy.int)
        else:
            self.TeffRange = numpy.array([0, 100000])
        if "loggRange" in keys:
            self.loggRange = numpy.array(self.config["loggRange"].split(),
                    dtype=numpy.float32)
        else:
            self.loggRange = numpy.array([2.5, 6.0])
        if "BfieldRange" in keys:
            self.BfieldRange = numpy.array(self.config["BfieldRange"].split(),
                    dtype=numpy.float32)
        else:
            self.BfieldRange = numpy.array([0.0, 20.0])
        self.resolvingPower = self.config["resolving_power"]
        if 'vsini' in keys:
            self.vsini = self.config['vsini']
        else:
            self.vsini = None
        self.Score = Score(directory=self.watchedDir)
        self.Score.selectMelodies(TeffRange=self.TeffRange, loggRange=self.loggRange,
                BfieldRange=self.BfieldRange, wlRange=self.wlRange)
        #, vsini=self.vsini,
        #        R=self.resolvingPower)
        #self.Score.setWlRange(self.wlStart, self.wlStop)

    def rehearse(self, vsini = None, R = None):
        if vsini == None:
            vsini = [self.vsini]
        if R == None:
            R = [self.resolvingPower]
        for v, r in zip(vsini, R):
            self.Score.rehearse(vsini=v, R=r)

    def perform(self, vsini=None, R = None, axes= None):

        if vsini == None:
            vsini = [self.vsini]
        if R == None:
            R = [self.resolvingPower]
        

        # LineStyles = different Vsini/R combinations
        # colors = different melodies (Teff/Log g/Bfield)
        keysignatures = []
        labels = []
        #linestyles = Lines.lineStyles.keys()[3:3+len(vsini)]
        linestyles = ['-', '--', '-.', ':']

        for v, r in zip(vsini, R):
            sp, l = self.Score.perform(vsini=v, R = r)
            keysignatures.append(sp)
            labels.append(l)

        if len(keysignatures) < 1:
            return

        colors = numpy.random.rand(len(keysignatures[0]), 3)
        for k, lb, ls in zip(keysignatures,labels, linestyles):
            for spectra, l, c in zip(k, lb, colors):
                for s in spectra:
                    line = axes.plot(s.wl, s.flux_I, ls=ls, color=c, label=l)

        #axes.figure.legend(lines, labels)

    def selectEnsemble(self, T=[], G=[], B=[]):
        self.Score.selectEnsemble(T=T, G=G, B=B)

    def inThePit(self):
        Teffs, loggs, Bfields = self.Score.getMelodyParams()
        return Teffs, loggs, Bfields

    def record(self, vsini=None, R=None):
        if vsini == None:
            vsini = [self.vsini]
        if R == None:
            R = [self.resolvingPower]

        for v, r, in zip(vsini, R):
            self.Score.record(vsini=v, R =r)
