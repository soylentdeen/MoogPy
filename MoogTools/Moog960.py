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
            HDUList.close()

class SyntheticPhrase( Phrase ):
    def __init__(self, rawData=[], interpolatedData=[], 
                 integratedData=[], convolvedData=[], diskInt = None):
        self.rawData = rawData
        if len(self.rawData) > 0:
            self.wlStart = rawData[0].header.get('WLSTART')
            self.wlStop = rawData[0].header.get('WLSTOP')
        elif len(interpolatedData) > 0:
            self.wlStart = interpolatedData[0].header.get('WLSTART')
            self.wlStop = interpolatedData[0].header.get('WLSTOP')
        elif len(integratedData) > 0:
            self.wlStart = integratedData[0].header.get('WLSTART')
            self.wlStop = integratedData[0].header.get('WLSTOP')
        elif len(convolvedData) > 0:
            self.wlStart = convolvedData[0].header.get('WLSTART')
            self.wlStop = convolvedData[0].header.get('WLSTOP')
        if diskInt == 'BEACHBALL':
            self.processedData = SpectralTools.BeachBall(parent=self,
                    interpolatedData=interpolatedData, integratedData=integratedData,
                    convolvedData=convolvedData)
        elif diskInt == 'DISKOBALL':
            self.processedData = SpectralTools.DiskoBall(parent=self,
                    interpolatedData=interpolatedData, integratedData=integratedData,
                    convolvedData=convolvedData)

    @classmethod
    def fromFile(self, hdr, data=None, filename=None, ext=None, diskInt=None,
                 sourceType="RAW"):
        if sourceType =="RAW":
            rawData = []
            rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
                filename=filename, ext=ext))
            interpolatedData = []
            integratedData = []
            convolvedData = []
        elif sourceType == "INTERPOLATED":
            rawData = []
            interpolatedData = []
            interpolatedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            integratedData = []
            convolvedData = []
        elif sourceType == "INTEGRATED":
            rawData = []
            interpolatedData = []
            integratedData = []
            integratedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            convolvedData = []
        elif sourceType == "CONVOLVED":
            rawData = []
            interpolatedData = []
            integratedData = []
            convolvedData = []
            convolvedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
        return self(rawData=rawData, interpolatedData=interpolatedData, 
                    integratedData=integratedData, convolvedData=convolvedData, diskInt=diskInt)
    
    def addSpectrum(self, hdr, data=None, filename=None, ext=None, sourceType="RAW"):
        if sourceType=="RAW":
            self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
        elif sourceType =="INTERPOLATED":
            self.processedData.interpolated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
        elif sourceType =="INTEGRATED":
            self.processedData.integrated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
        elif sourceType =="CONVOLVED":
            self.processedData.convolved.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))

    def tune(self, save=False, header=None):
        self.processedData.diskInt()
        if save:
            integratedFilename = self.rawData[0].filename[:-8]+'integrated.fits'
            self.saveIntegrated(filename = integratedFilename, header=header)

    def rehearse(self, vsini=0.0, R=0, observedWl=None):
        self.processedData.resample(vsini=vsini, R=R, observedWl=observedWl)

    def perform(self, vsini= 0.0, R = 0.0, observedWl = None, keySignature="CONVOLVED"):
        return self.processedData.yank(vsini=vsini, R = R, observedWl=observedWl,
                keySignature=keySignature)

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
            primary.header.set("SPECTRUM_CONTENTS", "RAW")
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)
            HDUList.close()


    def saveInterpolated(self, filename = None, header=None, primaryHeaderKWs={}):
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
            primary = pyfits.PrimaryHDU(header=header)
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            header.set("SPECTRUM_CONTENTS", "INTERPOLATED") 
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)
            HDUList.close()

    def saveIntegrated(self, filename = None, header=None, primaryHeaderKWs={}):
        HDUs = []
        for spectrum in self.processedData.integrated:
            hdr = spectrum.header.copy()
            spectrum.preserve(continuum=False)
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

            primary = pyfits.PrimaryHDU(header=header)
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            primary.header.set("SPECTRUM_CONTENTS", "INTEGRATED")
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)
            HDUList.close()

    def saveConvolved(self, vsini=None, R=None,filename = None, header=None, primaryHeaderKWs={}):
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
            primary = pyfits.PrimaryHDU(header=header)
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
            primary.header.set("SPECTRUM_CONTENTS", "CONVOLVED")
            HDUList.append(primary)
            for spectrum in HDUs:
                HDUList.append(spectrum)
            HDUList.update_extend()
            HDUList.verify(option='silentfix')
            HDUList.writeto(filename)
            HDUList.close()

class Melody( object ):
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        self.phrases = phrases
        self.header = header
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
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1]))

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
                self.phrases[i].saveConvolved(vsini=vsini, R=R, filename=filename, header=self.header)


class ObservedMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        super(ObservedMelody, self).__init__(phrases=phrases, filename=filename, label=label, header=header)

    @classmethod
    def fromFile(self, filename=None, label=None):
        info = pyfits.info(filename, output='')
        nPhrases = len(info)-1
        header = pyfits.getheader(filename, ext=0)
        phrases = []
        for i in range(nPhrases):
            hdr = pyfits.getheader(filename, ext=i+1)
            phrases.append(ObservedPhrase.fromFile(header=hdr, filename=filename, ext=i+1))

        return self(phrases=phrases, filename=filename, label=label, header=header)
        
    def loadData(self):
        for phrase in self.phrases:
            phrase.loadData()
            
    def rehearse(self, **kwargs):
        return

    def perform(self):
        spectra = []
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                spectra.append(self.phrases[i].listen())

        return spectra, self.label


class SyntheticMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        super(SyntheticMelody, self).__init__(phrases = phrases, filename=filename,
                label=label, header=header)
        self.loadMelody()

    def loadMelody(self):
        info = pyfits.info(self.filename, output='')
        self.nSpectra = len(info)-1
        self.header = pyfits.getheader(self.filename, ext=0, memmap=False)
        self.Teff = self.header.get("TEFF")
        self.logg = self.header.get("LOGG")
        self.B = self.header.get("BFIELD")
        self.contents = self.header.get("SPECTRUM_CONTENTS")
        if self.contents == None:
            self.contents = "RAW"
        self.rawLabel = "T = %dK log g = %.1f B = %.2f kG" % (self.Teff, self.logg, self.B)
        self.integratedLabels = []
        self.interpolatedLabels = []
        self.convolvedLabels = []

        self.phrases = []

        for i in range(self.nSpectra):
            added = False
            hdr = pyfits.getheader(self.filename, ext=i+1, memmap=False)
            for phrase in self.phrases:
                if phrase.owns(hdr):
                    phrase.addSpectrum(hdr, data=None, filename=self.filename,
                            ext=i+1, sourceType=self.contents)
                    added=True
                    break
            if not(added):
                self.phrases.append(SyntheticPhrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1, diskInt='BEACHBALL',
                    sourceType=self.contents))
                if self.contents == 'INTEGRATED':
                    self.integratedLabels.append("T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI')))
                if self.contents == 'CONVOLVED':
                    self.convolvedLabels.append("T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI'), hdr.get('RESOLVING_POWER')))

        self.nPhrases = len(self.phrases)

    def tune(self, save=False):
        for i in range(self.nPhrases):
            self.phrases[i].tune(save=save, header=self.header)

    def rehearse(self, vsini = 0.0, R = 0, observedWl = None):
        found = False
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].rehearse(vsini = vsini, R=R, observedWl=observedWl)
                found = True
        if found:
            self.convolvedLabels.append("T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" %
                    (self.Teff, self.logg, self.B, vsini, R))
            return self.convolvedLabels[-1]

    def perform(self, label="", keySignature="CONVOLVED"):
        params = {}
        if keySignature == "CONVOLVED":
            try:
                R = int(label.split('=')[5])
                vsini = float(label.split('=')[4].split()[0])
            except:
                R = 0
                vsini = 0.0
        elif keySignature == "INTEGRATED":
            R= 0
            try:
                vsini = float(label.split('=')[4].split()[0])
            except:
                vsini = 0.0
        spectra = []
        params["TEFF"] = self.Teff
        params["LOGG"] = self.logg
        params["BFIELD"] = self.B
        params["RESOLVING_POWER"] = R
        params["VSINI"] = vsini
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                spectra.append(self.phrases[i].perform(vsini=vsini, R=R, 
                    keySignature=keySignature))

        return spectra, params
    
class Score( object ):
    """
        This Score object contains many melodies.
    """
    def __init__(self, melodies = [], directory=None, observed=None, suffix='raw'):
        self.syntheticMelodies = melodies
        self.directory = directory
        self.observed = observed
        self.suffix = suffix
        self.loadMelodies()

    def loadMelodies(self):
        melodyFiles = glob.glob(self.directory+'*'+self.suffix+'.fits')
        for melody in melodyFiles:
            print("%s" % melody)
            self.syntheticMelodies.append(SyntheticMelody(filename=melody))

        if not(self.observed==None):
            self.Observed = ObservedMelody.fromFile(filename=self.observed, label='TWHydra')
            self.Observed.loadData()

    def getMelodyParams(self):
        raw_labels = []
        interpolated_labels = []
        integrated_labels = []
        convolved_labels = []
        observed_labels = []
        for melody in self.syntheticMelodies:
            raw_labels.append(melody.rawLabel)
            if len(melody.interpolatedLabels) > 0:
                for interpolated in melody.interpolatedLabels:
                    interpolated_labels.append(interpolated)
            if len(melody.integratedLabels) > 0:
                for integrated in melody.integratedLabels:
                    integrated_labels.append(integrated)
            if len(melody.convolvedLabels) > 0:
                for convolved in melody.convolvedLabels:
                    convolved_labels.append(convolved)
        if not(self.observed==None):
            observed_labels.append(self.Observed.label)

        return raw_labels, interpolated_labels, integrated_labels, convolved_labels, observed_labels

    def selectEnsemble(self, selectedLabels, keySignature='CONVOLVED'):
        if keySignature == 'RAW':
            for melody in self.syntheticMelodies:
                if melody.rawLabel in selectedLabels:
                    melody.muted=False
                else:
                    melody.muted=True
        elif keySignature == 'INTERPOLATED':
            for melody in self.syntheticMelodies:
                for interpolated in melody.interpolatedLabels:
                    if interpolated in selectedLabels:
                        melody.muted=False
                    else:
                        melody.muted=True
        elif keySignature == 'INTEGRATED':
            for melody in self.syntheticMelodies:
                for integrated in melody.integratedLabels:
                    if integrated in selectedLabels:
                        melody.muted=False
                    else:
                        melody.muted=True
        elif keySignature == 'CONVOLVED':
            for melody in self.syntheticMelodies:
                for convolved in melody.convolvedLabels:
                    if convolved in selectedLabels:
                        melody.muted=False
                    else:
                        melody.muted=True

    def addToEnsemble(self, selectedLabels):
        if keySignature == 'RAW':
            for melody in self.syntheticMelodies:
                if melody.rawLabel in selectedLabels:
                    melody.muted=False
        elif keySignature == 'INTERPOLATED':
            for melody in self.syntheticMelodies:
                for interpolated in melody.interpolatedLabels:
                    if interpolated in selectedLabels:
                        melody.muted=False
        elif keySignature == 'INTEGRATED':
            for melody in self.syntheticMelodies:
                for integrated in melody.integratedLabels:
                    if integrated in selectedLabels:
                        melody.muted=False
        elif keySignature == 'CONVOLVED':
            for melody in self.syntheticMelodies:
                for convolved in melody.convolvedLabels:
                    if convolved in selectedLabels:
                        melody.muted=False

    def selectMelodies(self, wlRange=[]):
        for melody in self.syntheticMelodies:
            #melody.muted=False
            #melody.inParameterRange(TeffRange=TeffRange,
            #    loggRange=loggRange, BfieldRange=BfieldRange)
            #if not(melody.muted):
            melody.selectPhrases(wlRange=wlRange)

        if not(self.observed==None):
            self.Observed.selectPhrases(wlRange=wlRange)

    def tune(self, save=False):
        for melody in self.syntheticMelodies:
            print melody
            melody.tune(save=save)


    def rehearse(self, vsini=0.0, R=0.0, binToObserved=True):
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
            if not(melody.muted):

                if binToObserved:
                    melody.rehearse(vsini=vsini, R=R, 
                        observedWl = self.compositeObserved.wl)
                else:
                    melody.rehearse(vsini=vsini, R=R, observedWl=None)

    def perform(self, selected=[], keySignature = "CONVOLVED"):
        spectra = []
        labels = []
        params = []
        for melody in self.syntheticMelodies:
            if not(melody.muted):
                if keySignature == "INTERPOLATED":
                    for i in range(len(melody.interpolatedLabels)):
                        if melody.interpolatedLabels[i] in selected:
                            sp, p = melody.perform(melody.interpolatedLabels[i],
                                    keySignature=keySignature)
                            spectra.append(sp)
                            params.append(p)
                            labels.append(melody.interpolatedLabels[i])

                if keySignature == "INTEGRATED":
                    for i in range(len(melody.integratedLabels)):
                        if melody.integratedLabels[i] in selected:
                            sp, p = melody.perform(melody.integratedLabels[i],
                                    keySignature=keySignature)
                            spectra.append(sp)
                            params.append(p)
                            labels.append(melody.integratedLabels[i])

                if keySignature == "CONVOLVED":
                    for i in range(len(melody.convolvedLabels)):
                        if melody.convolvedLabels[i] in selected:
                            sp, p = melody.perform(melody.convolvedLabels[i],
                                    keySignature=keySignature)
                            spectra.append(sp)
                            params.append(p)
                            labels.append(melody.convolvedLabels[i])

        return spectra, params, labels


    def listen(self):   # gets the observed spectrum
        spectra, label = self.Observed.perform()
        compositeSpectrum = None
        for sp in spectra:
            compositeSpectrum = SpectralTools.mergeSpectra(first=compositeSpectrum,
                    second=sp)
        
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
