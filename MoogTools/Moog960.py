import pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy
import os
import matplotlib.lines as Lines

class Phrase( object ):
    def __init__(self, wlStart, wlStop):
        self.wlStart = wlStart
        self.wlStop = wlStop
        return

    def owns(self, hdr):
        if ((self.wlStart == hdr.get('WLSTART')) & (self.wlStop == hdr.get("WLSTOP"))):
            return True
        return False

    def inWlRange(self, wlStart, wlStop):
        print self.wlStart, wlStart
        print self.wlStop, wlStop
        return ((self.wlStart < wlStop) & (self.wlStop > wlStart))

class ObservedPhrase( Phrase ):
    def __init__(self, observedData = None):
        observedData = SpectralTools.ObservedSpectrum(observed=observedData)
        wlStart = observedData.observed.header.get('WLSTART')
        wlStop = observedData.observed.header.get('WLSTOP')
        super(ObservedPhrase, self).__init__(wlStart, wlStop)

        self.observedData = observedData

    @classmethod
    def fromFile(self, header=None, data=None, filename=None, ext=None):
        observed = SpectralTools.Spectrum.from_file(header= header, data=data,
                filename=filename, ext=ext)
        return self(observedData=observed)
        
    def loadData(self):
        self.observedData.observed.loadData()
        
    def listen(self):
        return self.observedData.observed
        
    def record(self, filename = None):
        self.save(filename=filename)

    def save(self, filename = None):
        HDUs = []
        spectrum = self.observedData.observed
        hdr = spectrum.header.copy()
        SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns, header=hdr)
        SpectrumHDU.name = "%.4fA - %.4fA OBSERVED" % (hdr.get('wlStart'), 
                hdr.get('wlStop'))
        HDUs.append(SpectrumHDU)

        if filename == None:
            return HDUs

        if os.path.exists(filename):    #file exists!  get ready to append!
            print filename
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
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)

class SyntheticPhrase( Phrase ):
    def __init__(self, rawData=None, diskInt = None):
        self.rawData = rawData
        self.wlStart = rawData[0].header.get('WLSTART')
        self.wlStop = rawData[0].header.get('WLSTOP')
        if diskInt == 'BEACHBALL':
            self.processedData = SpectralTools.BeachBall(parent=self)
        elif diskInt == 'DISKOBALL':
            self.processedData = SpectralTools.DiskoBall(parent=self)

    @classmethod
    def fromFile(self, hdr, data=None, filename=None, ext=None, diskInt=None):
        rawData = []
        rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
            filename=filename, ext=ext))
        return self(rawData=rawData, diskInt=diskInt)
    
    def addRawSpectrum(self, hdr, data=None, filename=None, ext=None):
        self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
            filename=filename, ext=ext))

    def rehearse(self, vsini=0.0, R=0, observedWl=None):
        self.processedData.resample(vsini=vsini, R=R, observedWl=observedWl)

    def perform(self, vsini= 0.0, R = 0.0, observedWl = None):
        return self.processedData.yank(vsini=vsini, R = R, observedWl=observedWl)

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

    def saveConvolved(self, vsini=None, R=None,filename = None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.processedData.convolved:
            hdr = spectrum.header.copy()
            if (not(vsini==None) and not(R==None)):
                if (vsini == hdr.get('VSINI')) and (R == hdr.get('RESOLVING_POWER')):
                    spectrum.preserve(continuum=False)
                    SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                        header=hdr)
                    SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f R=%.1f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'), 
                        hdr.get('RESOLVING_POWER'))
            
                    HDUs.append(SpectrumHDU)
            else:
                spectrum.preserve(continuum=False)
                SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
                SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f R=%.1f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'), 
                    hdr.get('RESOLVING_POWER'))
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
    def __init__(self, phrases = [], filename=None, label=None):
        self.phrases = phrases
        self.nPhrases = len(self.phrases)
        self.filename = filename
        self.selectedPhrases = [False for i in range(self.nPhrases)]
        self.muted = True
        self.label = label
        
    def addPhrase(self, phrases):
        for phrase in phrases:
            self.phrases.append(phrase)
        self.nPhrases = len(self.phrases)

    def selectPhrases(self, wlRange=[], selectAll=False):
        self.selectedPhrases = []
        print "Junk"
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                print "Blah"
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1]))
        print self.selectedPhrases

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
    

    def record(self, label="", basename = None):
        R = int(label.split('=')[5])
        vsini = float(label.split('=')[4].split()[0])
        print("Recording record \'%s\' to disk" % label)
        if basename == None:
            filename = self.filename[:-4]+'_saved.fits'
        else:
            filename = basename+"_T%d_G%.2f_B%.2f_R%d_V%.2f.fits" % (
                self.Teff, self.logg, self.B, R, vsini)
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].saveConvolved(vsini=vsini, R=R, filename=filename)


class ObservedMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None):
        super(ObservedMelody, self).__init__(phrases=phrases, filename=filename, label=label)

    @classmethod
    def fromFile(self, filename=None, label=None):
        info = pyfits.info(filename, output='')
        nPhrases = len(info)-1
        header = pyfits.getheader(filename, ext=0)
        phrases = []
        for i in range(nPhrases):
            hdr = pyfits.getheader(filename, ext=i+1)
            phrases.append(ObservedPhrase.fromFile(header=hdr, filename=filename, ext=i+1))

        return self(phrases=phrases, filename=filename, label=label)
        
    def loadData(self):
        for phrase in self.phrases:
            phrase.loadData()
            
    def rehearse(self, **kwargs):
        return

    def perform(self):
        spectra = []
        for i in range(self.nPhrases):
            print self.selectedPhrases[i]
            if self.selectedPhrases[i]:
                spectra.append(self.phrases[i].listen())

        print spectra
        print self.label
        return spectra, self.label


class SyntheticMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None):
        super(SyntheticMelody, self).__init__(phrases = phrases, filename=filename,
                label=label)
        self.loadMelody()

    def loadMelody(self):
        info = pyfits.info(self.filename, output='')
        self.nSpectra = len(info)-1
        self.header = pyfits.getheader(self.filename, ext=0)
        self.Teff = self.header.get("TEFF")
        self.logg = self.header.get("LOGG")
        self.B = self.header.get("BFIELD")
        self.label = ["T = %dK log g = %.1f B = %.2f kG" % (self.Teff, self.logg, self.B)]

        self.phrases = []

        for i in range(self.nSpectra):
            added = False
            hdr = pyfits.getheader(self.filename, ext=i+1)
            for phrase in self.phrases:
                if phrase.owns(hdr):
                    phrase.addRawSpectrum(hdr, data=None, filename=self.filename,
                            ext=i+1)
                    added=True
                    break
            if not(added):
                self.phrases.append(SyntheticPhrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1, diskInt='BEACHBALL'))

        self.nPhrases = len(self.phrases)

    def rehearse(self, vsini = 0.0, R = 0, observedWl = None):
        found = False
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].rehearse(vsini = vsini, R=R, observedWl=observedWl)
                found = True
        if found:
            self.label.append("T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % (self.Teff, self.logg, self.B, vsini, R))

    def perform(self, label):
        R = int(label.split('=')[5])
        vsini = float(label.split('=')[4].split()[0])
        spectra = []
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                spectra.append(self.phrases[i].perform(vsini=vsini, R=R))

        return spectra
    
class Score( object ):
    """
        This Score object contains many melodies.
    """
    def __init__(self, melodies = [], directory=None, observed=None):
        self.syntheticMelodies = melodies
        self.directory = directory
        self.observed = observed
        self.loadMelodies()

    def loadMelodies(self):
        melodyFiles = glob.glob(self.directory+'*raw.fits')
        for melody in melodyFiles:
            print("%s" % melody)
            m = SyntheticMelody(filename=melody)
            self.syntheticMelodies.append(m)

        if not(self.observed==None):
            self.Observed = ObservedMelody.fromFile(filename=self.observed, label='TWHydra')
            self.Observed.loadData()

    def getMelodyParams(self):
        Teff = []
        logg = []
        B = []
        raw_labels = []
        processed_labels = []
        for melody in self.syntheticMelodies:
            Teff.append(melody.Teff)
            logg.append(melody.logg)
            B.append(melody.B)
            raw_labels.append(melody.label[0])
            print "length of melody %d" % len(melody.label)
            print melody.label
            if len(melody.label) > 1:
                for convolved in melody.label[1:]:
                    processed_labels.append(convolved)

        return raw_labels, processed_labels

    def selectEnsemble(self, selectedLabels):
        for melody in self.syntheticMelodies:
            if melody.label[0] in selectedLabels:
                melody.muted=False
                print("Loud: %d, %.1f, %.1f" % (melody.Teff, melody.logg, melody.B))
            else:
                melody.muted=True
                print("Mute: %d, %.1f, %.1f" % (melody.Teff, melody.logg, melody.B))

    def addToEnsemble(self, selectedLabels):
        print "Attempting to add labels"
        print selectedLabels
        for melody in self.syntheticMelodies:
            if melody.label[0] in selectedLabels:
                melody.muted=False
                print("Loud: %d, %.1f, %.1f" % (melody.Teff, melody.logg, melody.B))


    def selectMelodies(self, wlRange=[]):
        for melody in self.syntheticMelodies:
            #melody.muted=False
            #melody.inParameterRange(TeffRange=TeffRange,
            #    loggRange=loggRange, BfieldRange=BfieldRange)
            #if not(melody.muted):
            melody.selectPhrases(wlRange=wlRange)

        self.Observed.selectPhrases(wlRange=wlRange)

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
        for melody in self.syntheticMelodies:
            print("Score::Melody::Rehearse: %s" %melody.muted)
            if not(melody.muted):
                melody.rehearse(vsini=vsini, R=R, 
                        observedWl = self.compositeObserved.wl)

    def perform(self, selected=[]):
        spectra = []
        labels = []
        for melody in self.syntheticMelodies:
            if not(melody.muted):
                print("Aha! We have a non-muted melody!")
                for i in range(len(melody.label)-1):
                    if melody.label[i+1] in selected:
                        print("Trying a performance!")
                        spectra.append(melody.perform(melody.label[i+1]))
                        labels.append(melody.label[i+1])

        return spectra, labels


    def listen(self):   # gets the observed spectrum
        spectra, label = self.Observed.perform()
        compositeSpectrum = None
        for sp in spectra:
            compositeSpectrum = SpectralTools.mergeSpectra(first=sp,
                    second=compositeSpectrum)
        
        self.compositeObserved = compositeSpectrum
        return compositeSpectrum, label

    def record(self, selected=[], basename=''):
        for melody in self.syntheticMelodies:
            if not(melody.muted):
                for i in range(len(melody.label)-1):
                    if melody.label[i+1] in selected:
                        melody.record(melody.label[i+1], basename=basename)

class Moog960( object ):
    def __init__(self, configFile):
        """
        Moog960 - A Mixer for synthetic/observed spectra
        
        Usage:
                mixer = Moog960(configFile)
        
        Member Functions:
            Moog960.applyConfigFile() - private
                  - applies the configuration file
            Moog960.rehearse(vsini=None, R=None)
                  - Ensures that the supplied rotational broadening and
                  resolvining power has been applied to all synthetic meodies.
                  If not, the synthesizer will trigger these calcuation
            Moog960.perform(vsini=None, R=None, plotaxes=None, diffaxes=None)
                  - For the 
                    
            
        """
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

        if 'observed' in keys:
            self.observed = self.config['observed']
            
        if 'wlShift' in keys:
            self.wlShift = self.config['wlShift']
        else:
            self.wlShift = 0.0

        self.Score = Score(directory=self.watchedDir, observed=self.observed)
        #self.Score.selectMelodies(TeffRange=self.TeffRange, loggRange=self.loggRange,
        #        BfieldRange=self.BfieldRange, wlRange=self.wlRange)
        self.Score.selectMelodies(wlRange=self.wlRange)
        mainTheme, obs_label = self.Score.listen()
        mainTheme.wl += self.wlShift
        self.mainTheme = mainTheme
        self.observed_label = obs_label

    def rehearse(self, vsini = None, R = None):
        if vsini == None:
            vsini = [self.vsini]   # Should this permit a list of vsinis?
        if R == None:
            R = [self.resolvingPower]
        for v, r in zip(vsini, R):
            self.Score.rehearse(vsini=v, R=r)

    
    #def perform(self, vsini=None, R = None, plotaxes= None, diffaxes = None):
    def perform(self, selected=[], plotaxes=None, diffaxes=None):

        labels = []
        performances = []
        linestyles = ['-', '--', '-.', ':']

        performances, labels = self.Score.perform(selected=selected)
        
        colors = numpy.random.rand(len(performances), 3)
        print performances
        for spectrum, l, c in zip(performances, labels, colors):
            for s in spectrum:
                if not(plotaxes==None):
                    line = plotaxes.plot(s.wl, s.flux_I, color=c, label=l)
                if not(diffaxes==None):
                    difference = s - self.mainTheme
                    line = diffaxes.plot(difference.wl, difference.flux_I, color=c, label=l)
                
        #observed
        plotaxes.plot(self.mainTheme.wl, self.mainTheme.flux_I, color='k')
        
        plotaxes.set_xbound(lower=self.wlRange[0], upper=self.wlRange[1])
        if not(diffaxes==None):
            diffaxes.set_xbound(lower=self.wlRange[0], upper=self.wlRange[1])
        #axes.figure.legend(lines, labels)

    def selectEnsemble(self, selected):
        self.Score.selectEnsemble(selected)
        
    def addToEnsemble(self, selected):
        self.Score.addToEnsemble(selected)

    def inThePit(self):
        melodies, processed = self.Score.getMelodyParams()
        return melodies, processed

    def record(self, selected=[], basename=""):

        self.Score.record(selected=selected, basename=basename)
