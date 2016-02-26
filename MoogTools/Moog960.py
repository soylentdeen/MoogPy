import pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy
import os
import matplotlib.lines as Lines
import time

class Label( object ):
    def __init__(self, parameters):
        self.parameters = parameters

    def copy(self):
        return Label(self.parameters)

    def __eq__(self, other):
        for key in self.parameters.keys():
            if self.parameters[key] != other.parameters[key]:
                return False
        return True
        #return ((self.parameters["LABEL"] == other.parameters["LABEL"]) &
        #        (self.parameters["WLSTART"] == other.parameters["WLSTART"]) &
        #        (self.parameters["WLSTOP"] == other.parameters["WLSTOP"]))

    def __cmp__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart().__cmp__(other.getWlStart())

    def getWlStart(self):
        return self.parameters["WLSTART"]

class Moog960Error( Exception ):
    def __init__(self, value, errmsg):
        self.value = value
        self.message = {}
        self.message[0] = "Failure retriving Raw Data!! %s" % errmsg
        self.message[1] = "Failure loading Processed Data!! %s" % errmsg

    def __str__(self):
        return repr(self.message[self.value])


class Phrase( object ):
    """
    Moog960::Phrase
    
    A phrase contains the output of an individual synthesis (unique Teff, log g, and B field)
        over a single wavelength region, both raw and processed data.
    
    """
    def __init__(self, wlStart, wlStop):
        self.wlStart = wlStart
        self.wlStop = wlStop
        return

    def owns(self, hdr):
        """
        Moog960::Phrase:owns(self, pyfitsHeader):
        
        This routine checks to see if a fits Header belongs to this Phrase
        
        Input:
            hdr - pyfits Header object
        
        Returns:
            True : If the phrase's wlStart and wlStop are equal to the wlStart/Stop in the header
            False: If not.
        """
        if ((self.wlStart == hdr.get('WLSTART')) & (self.wlStop == hdr.get("WLSTOP"))):
            return True
        return False

    def inWlRange(self, wlStart, wlStop):
        """
        Moog960::Phrase::inWlRange(self, wlStart, wlStop)
        
        This routine checks if the phrase has any overlap with the given wavelength region
        
        Input:
            wlStart - [float] - Start of wavelength region
            wlStop - [float] - End of wavelength region
            
        Returns:
            True : if the Phrase overlaps with the wavelength region
            False: if not.
        """
        return ((self.wlStart < wlStop) & (self.wlStop > wlStart))

class ObservedPhrase( Phrase ):
    def __init__(self, observedData = [], observedLabels = []):
        wlStart = observedData[0].header.get('WLSTART')
        wlStop = observedData[0].header.get('WLSTOP')
        super(ObservedPhrase, self).__init__(wlStart, wlStop)

        self.observedData = observedData
        self.observedLabels = observedLabels

    @classmethod
    def fromFile(self, header=None, data=None, filename=None, ext=None):
        observed = SpectralTools.Spectrum.from_file(header= header, data=data,
                filename=filename, ext=ext)
        parameters = {}
        parameters["WLSTART"] = observed.header.get('WLSTART')
        parameters["WLSTOP"] = observed.header.get('WLSTOP')
        return self(observedData=[observed], observedLabels = [Label(parameters)])
        
    def loadData(self):
        for observed in self.observedData:
            observed.loadData()
        
    def listen(self):
        return self.observedData, self.observedLabels
        
    def record(self, filename = None):
        self.save(filename=filename)

    def save(self, filename = None):
        HDUs = []
        for spectrum in self.observedData:
            hdr = spectrum.header.copy()
            spectrum.preserve(continuum=False, V=False)
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
                    print spectrum.name
                    print HDUList
                    HDUList.pop(HDUList.index_of(spectrum.name))
                except:
                    pass
                HDUList.append(spectrum)
                print len(HDUList)
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
                 integratedData=[], convolvedData=[], 
                 rawLabels=[], interpolatedLabels=[],
                 integratedLabels=[], convolvedLabels=[], diskInt = None):
        """
        Moog960::SyntheticPhrase::init(rawData=[], interpolatedData=[], integratedData=[]
                                       convolvedData=[], diskInt=None):
        Creates a Synthetic phrase from synthesized MoogStokes Data
        
        Input:
            rawData = list of raw data Spectrum objects containing the raw data
            interpolatedData = list of interpolated Spectrum objects
            integratedData = list of integrated Spectrum objects
            convolvedData = list of convolved Spectrum objects
        """
        self.rawData = rawData
        self.rawLabels = rawLabels
        self.interpolatedData = interpolatedData
        self.interpolatedLabels = interpolatedLabels
        self.integratedData = integratedData
        self.integratedLabels = integratedLabels
        self.convolvedData = convolvedData
        self.convolvedLabels = convolvedLabels
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
            self.processedData = SpectralTools.BeachBall(parent=self)
        elif diskInt == 'DISKOBALL':
            self.processedData = SpectralTools.DiskoBall(parent=self)

    @classmethod
    def fromFile(self, hdr, data=None, filename=None, ext=None, diskInt=None,
                 sourceType="RAW", parameters=None):
        """
        Creates a phrase from a data file    
        """
        if sourceType =="RAW":
            rawData = []
            rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
                filename=filename, ext=ext))
            interpolatedData = []
            integratedData = []
            convolvedData = []
            
            rawLabels = []
            rawLabels.append(Label(parameters))
            interpolatedLabels = []
            integratedLabels = []
            convolvedLabels = []
        elif sourceType == "INTERPOLATED":
            rawData = []
            interpolatedData = []
            interpolatedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            integratedData = []
            convolvedData = []

            rawLabels = []
            interpolatedLabels = []
            interpolatedLabels.append(Label(parameters))
            integratedLabels = []
            convolvedLabels = []
        elif sourceType == "INTEGRATED":
            rawData = []
            interpolatedData = []
            integratedData = []
            integratedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            convolvedData = []
            
            rawLabels = []
            interpolatedLabels = []
            integratedLabels = []
            integratedLabels.append(Label(parameters))
            convolvedLabels = []
        elif sourceType == "CONVOLVED":
            rawData = []
            interpolatedData = []
            integratedData = []
            convolvedData = []
            convolvedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
                
            rawLabels = []
            interpolatedLabels = []
            integratedLabels = []
            convolvedLabels = []
            convolvedLabels.append(Label(parameters))

        return self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                rawLabels=rawLabels, interpolatedLabels=interpolatedLabels,
                integratedLabels=integratedLabels, convolvedLabels=convolvedLabels, diskInt=diskInt)
    
    def addSpectrum(self, hdr, data=None, filename=None, ext=None, sourceType="RAW", parameters=None):
        """
        Moog960::SyntheticPhrase::addSpectrum(hdr, data=None, filename=None, ext=None, sourceType="RAW")
        
        Adds a spectrum to the phrase, depending on its source (RAW, INTERPOLATED, INTEGRATED, CONVOLVED)
        """
        if sourceType=="RAW":
            self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.rawLabels.append(Label(parameters))
        elif sourceType =="INTERPOLATED":
            self.processedData.interpolated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.interpolatedLabels.append(Label(parameters))
        elif sourceType =="INTEGRATED":
            self.processedData.integrated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.integratedLabels.append(Label(parameters))
        elif sourceType =="CONVOLVED":
            self.processedData.convolved.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.convolvedLabels.append(Label(parameters))

    def tune(self, vsini=0.0, save=False, header=None):
        """
        Moog960::SyntheticPhrase::tune(vsini=0.0, save=False, header=None)
        
        tune() is the first step in processing raw synthetic spectra.  It performs the
            disk integration and convolution with a rotational broadening kernel.
        
        Input:
            vsini = [float] rotational broadening (km/s)
            save = True/False - whether or not to save the disk integrated spectrum
            header = optional pyfits header to be saved with the integrated spectrum
        """
        created = self.processedData.diskInt(vsini=vsini)
        if created:
            parameters = self.rawLabels[0].parameters.copy()
            parameters["VSINI"] = vsini
            self.integratedLabels.append(Label(parameters))
        if save:
            integratedFilename = self.rawData[0].filename[:-8]+'integrated.fits'
            self.saveIntegrated(filename = integratedFilename, header=header)
            
        # TODO: if a new Integrated Spectrum is created in processedData, update integratedLabels list

    def rehearse(self, vsini=0.0, R=0.0, observedWl=None):
        """
        Moog960::SyntheticPhrase::rehearse(vsini=0.0, R=0, observedWl=None)
        
        rehearse() is the final step in producing synthetic spectra suitable for 
            comparison with observed spectra
            
        Input:
            vsini = (float) rotational broadening (km/s)
            R = (float) resolving power of convolved data
            observedWl =  numpy.array([float]) - observed wavelength
            
        Returns:
            Nothing
        """
        created = self.processedData.resample(vsini=vsini, R=R, observedWl=observedWl)
        if 'INTEGRATED' in created:
            parameters = self.rawLabels[0].parameters.copy()
            parameters["VSINI"] = vsini
            parameters["WLSTART"] = self.wlStart
            parameters["WLSTOP"] = self.wlStop
            parameters["LABEL"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini)
            self.integratedLabels.append(Label(parameters))
        if 'CONVOLVED' in created:
            parameters = self.rawLabels[0].parameters.copy()
            parameters["VSINI"] = vsini
            parameters["R"] = R
            parameters["WLSTART"] = self.wlStart
            parameters["WLSTOP"] = self.wlStop
            parameters["LABEL"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini, R)

            self.convolvedLabels.append(Label(parameters))

    #def perform(self, vsini= 0.0, R = 0.0, observedWl = None, keySignature="CONVOLVED"):
    def perform(self, label=None, keySignature="CONVOLVED"):
        """
        Moog960::SyntheticPhrase::perform(vsini=0.0, R=0.0, observedWl=None, keySignature="CONVOLVED")
        
        perform() returns a previously calculated spectrum for further processing or display
        
        Input:
            vsini = (float) rotational broadening (km/s)
            R = (float) resolving power of convolved data
            observedWl =  numpy.array([float]) - observed wavelength
            keySignature = "RAW", "INTERPOLATED", "INTEGRATED", or "CONVOLVED"
            
        Output:
            Spectrum corresponding to the desired spectrum
        """
        #TODO instead of yank, examine SyntheticPhrase labels
        
        if keySignature == "CONVOLVED":
            for l, convolved in zip(self.convolvedLabels, self.convolvedData):
                if l==label:
                    return convolved, label
        
        raise Moog960Error (1, "The requested spectrum does not exist in the selected Phrase!")
        
        #return self.processedData.yank(vsini=vsini, R = R, observedWl=observedWl,
        #        keySignature=keySignature)

    def saveRaw(self, filename = None, primaryHeaderKWs={}):
        """
        Moog960::SyntheticPhrase::saveRaw(filename = None, primaryHeaderKWs={})
        
        Saves raw data to a file
        
        Input:
            filename = name of file to which data will be saved
            primaryHeaderKWs = dictionary containing fits header keywords and values
            
        Output:
            Nothing
        """
        HDUs = []
        for spectrum in self.rawData:
            hdr = spectrum.header.copy()
            SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                    header=hdr)
            SpectrumHDU.name = "%.4fA - %.4fA PHI=%.3f MU=%.3f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get("PHI_ANGLE"), hdr.get("MU"))
            
            #TODO: append relevant data from the SyntheticPhrase.rawLabel

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
        """
        Moog960::SyntheticPhrase::saveInterpolated(filename = None, primaryHeaderKWs={})
        
        Saves interpolated data to a file
        
        Input:
            filename = name of file to which data will be saved
            primaryHeaderKWs = dictionary containing fits header keywords and values
            
        Output:
            Nothing
        """
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
        """
        Moog960::SyntheticPhrase::saveIntegrated(filename = None, primaryHeaderKWs={})
        
        Saves integrated data to a file
        
        Input:
            filename = name of file to which data will be saved
            primaryHeaderKWs = dictionary containing fits header keywords and values
            
        Output:
            Nothing
        """
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

    def saveConvolved(self, vsini=None, R=None,filename = None, header=None, 
            wlStart=None, wlStop=None, label=None, primaryHeaderKWs={}):
        """
        Moog960::SyntheticPhrase::saveConvolved(filename = None, primaryHeaderKWs={}
                vsini=None, R=None, header=None, wlStart=None, wlStop=None)
        
        Saves convolved data to a file
        
        Input:
            filename = name of file to which data will be saved
            primaryHeaderKWs = dictionary containing fits header keywords and values
            
        Output:
            Nothing
        """

        HDUs = []
        if label != None:
            try:
                index = self.convolvedLabels.index(label)
                spectrum = self.processedData.convolved[index]
                hdr = spectrum.header.copy()
                spectrum.preserve(continuum=False)
                SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                        header = hdr)
                SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f R=%.1f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'), 
                                            hdr.get('RESOLVING_POWER'))
                HDUs.append(SpectrumHDU)
            except:
                raise Moog960Error (1, "The requested spectrum does not exist in the selected Phrase!")
        else:
            for spectrum in self.processedData.convolved:
                hdr = spectrum.header.copy()
                # TODO : implement Label stuff
                if (not(vsini==None) and not(R==None) and not(wlStart==None) and
                        not(wlStop==None)):
                    if ((vsini == hdr.get('VSINI')) and (R == hdr.get('RESOLVING_POWER')) 
                            and (wlStart == hdr.get('WLSTART')) and (wlStop == hdr.get('WLSTOP'))):
                        spectrum.preserve(continuum=False)
                        SpectrumHDU = pyfits.BinTableHDU.from_columns(spectrum.columns,
                            header=hdr)
                        SpectrumHDU.name = "%.4fA - %.4fA VSINI=%.3f R=%.1f" % (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get('VSINI'), 
                            hdr.get('RESOLVING_POWER'))
                
                        HDUs.append(SpectrumHDU)

                elif (not(vsini==None) and not(R==None)):
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
    """
    Moog960::Melody
    
    A melody object can contain many phrases which correspond to the same physical 
        parameters (Teff, log g, and B field), but different wavelength regions.
    
    Member variables:
        phrases = list [] of phrases
        header = Primary fits header for the melody (Teff, logg, B)
        nPhrases = number of phrases in the melody
        filename = name of the .fits file associated with the melody
        selectedPhrases = list [] of Booleans signifying whether or not a phrase is selected - obsolete?
        muted = Boolean signifiying whether the melody as a whole is selected - obsolete?
        label = text label describing the melody (Teff, logg, B) - obsolete?
    
    """
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        self.phrases = phrases
        self.header = header
        self.nPhrases = len(self.phrases)
        self.filename = filename
        self.selectedPhrases = [False for i in range(self.nPhrases)]
        self.muted = True
        self.label = label
        
    def addPhrase(self, phrases):    # WTF?  This doesn't look right
        for phrase in phrases:
            self.phrases.append(phrase)
        self.nPhrases = len(self.phrases)

    def selectPhrases(self, wlRange=[], selectAll=False):
        """
        Moog960::Melody::selectPhrases(wlRange=[], selectAll=False):
        
        selectPhrases() goes through the phrases in a melody and selects
            those phrases which overlap with the set wavelength range
            
        Input:
            wlRange = [wlStart, wlStart]
            
        Output:
            Nothing
        """
        self.selectedPhrases = []
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1]))

    def record(self, labels=None, basename = None):
        """
        Moog960::Melody::record(self, labels=[], basename=None)
        
        record() saves the selected phrases of a melody to disk
        """
        
        filename = self.filename[:-4]+'_saved.fits'
        if labels != None:
            for label in labels:
                R = label.parameters["R"]
                vsini = label.parameters["VSINI"]
                if basename != None:
                    filename = basename+"_T%d_G%.2f_B%.2f_R%d_V%.2f.fits" % (
                        self.Teff, self.logg, self.B, R, vsini)
                print("Recording record \'%s\' to disk" % label.parameters["LABEL"])
                for i in range(self.nPhrases):
                    if self.selectedPhrases[i]:
                       self.phrases[i].saveConvolved(vsini=vsini, R=R,
                               filename=filename, header=self.header, 
                               wlStart=label.parameters["WLSTART"], 
                               wlStop = label.parameters["WLSTOP"])
        else:
            for i in range(self.nPhrases):
                if self.selectedPhrases[i]:
                    self.phrases[i].saveConvolved()

                           
            #TODO make Labels actually do what they are supposed to.


class ObservedMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        super(ObservedMelody, self).__init__(phrases=phrases, filename=filename, 
                label=label, header=header)

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

    def selectPhrases(self, wlRange=[], selectAll=False):
        self.selectedPhrases = []
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1]))

    def perform(self):
        spectra = []
        labels = []
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                sp, l = self.phrases[i].listen()
                spectra.append(sp)
                labels.append(l)

        return spectra, labels

    def record(self, labels=[], basename = None):
        """
        Moog960::Melody::record(self, labels=[], basename=None)
        
        record() saves the selected phrases of a melody to disk
        """
        for label in labels:
            print("Recording record \'%s\' to disk" % label.parameters["LABEL"])
            if basename == None:
                filename = self.filename[:-4]+'_saved.fits'
            else:
                filename = basename+".fits"
            for i in range(self.nPhrases):
                if (self.selectedPhrases[i] and (label in
                    self.phrases[i].observedLabels)):
                   print label.parameters
                   self.phrases[i].save(filename=filename)
                           

class SyntheticMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None, header=None):
        super(SyntheticMelody, self).__init__(phrases = phrases, filename=filename,
                label=label, header=header)
        if len(self.phrases) == 0:
            self.loadMelody()
        else:
            self.Teff = self.header.get("TEFF")
            self.logg = self.header.get("LOGG")
            self.B = self.header.get("BFIELD")
            self.contents = self.header.get("SPECTRUM_CONTENTS")

    def loadMelody(self):
        info = pyfits.info(self.filename, output='')
        nSpectra = len(info)-1
        self.header = pyfits.getheader(self.filename, ext=0, memmap=False)
        self.phrases = []
        self.Teff = self.header.get("TEFF")
        self.logg = self.header.get("LOGG")
        self.B = self.header.get("BFIELD")
        self.contents = self.header.get("SPECTRUM_CONTENTS")
        if self.contents == None:
            self.contents = "RAW"
        parameters = {}
        parameters["LABEL"] = "T = %dK log g = %.1f B = %.2f kG" % (self.Teff, self.logg, self.B)
        parameters["TEFF"] = self.Teff
        parameters["LOGG"] = self.logg
        parameters["BFIELD"] = self.B

        for i in range(nSpectra):
            added = False
            hdr = pyfits.getheader(self.filename, ext=i+1, memmap=False)
            for phrase in self.phrases:
                if phrase.owns(hdr):
                    phrase.addSpectrum(hdr, data=None, filename=self.filename,
                        ext=i+1, sourceType=self.contents, parameters=parameters)
                    added=True
                    break
            if not(added):
                parameters = {}
                parameters["TEFF"] = self.Teff
                parameters["LOGG"] = self.logg
                parameters["BFIELD"] = self.B
                if self.contents == 'INTEGRATED':
                    label = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI'))
                    parameters["LABEL"] = label
                    parameters["VSINI"] = hdr.get('VSINI')
                    parameters["WLSTART"] = hdr.get('WLSTART')
                    parameters["WLSTOP"] = hdr.get('WLSTOP')
                    parameters["SELECTED"] = False
                if self.contents == 'CONVOLVED':
                    label = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI'), hdr.get('RESOLVING_POWER'))
                    parameters["LABEL"] = label
                    parameters["VSINI"] = hdr.get('VSINI')
                    parameters["R"] = hdr.get('RESOLVING_POWER')
                    parameters["WLSTART"] = hdr.get('WLSTART')
                    parameters["WLSTOP"] = hdr.get('WLSTOP')
                    parameters["SELECTED"] = False
                self.phrases.append(SyntheticPhrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1, diskInt='BEACHBALL',
                    sourceType=self.contents, parameters=parameters))


        self.nPhrases = len(self.phrases)

    def tune(self, save=False):
        for i in range(self.nPhrases):
            self.phrases[i].tune(save=save, header=self.header)

    def selectPhrases(self, wlRange=[], selectAll=False, keySignature="CONVOLVED"):
        self.selectedPhrases = []
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1]))

    def rehearse(self, vsini = 0.0, R = 0, observedWl = None, returnLabels =False):
        convolved = []
        for i in range(self.nPhrases):
            if self.selectedPhrases[i]:
                self.phrases[i].rehearse(vsini = vsini, R=R, observedWl=observedWl)
                convolved.append(self.phrases[i].convolvedLabels[-1])
        if returnLabels:
            return convolved

    def perform(self, label=None, keySignature="CONVOLVED"):
        """
        if keySignature == "CONVOLVED":
            try:
                R = label.parameters["R"]
                vsini = label.parameters["VSINI"]
            except:
                R = 0
                vsini = 0.0
        elif keySignature == "INTEGRATED":
            R= 0
            try:
                vsini = label.parameters["VSINI"]
            except:
                vsini = 0.0
        """
        for i in range(self.nPhrases):
            if (self.selectedPhrases[i] & (label.parameters["WLSTART"] == self.phrases[i].wlStart) &
                    (label.parameters["WLSTOP"] == self.phrases[i].wlStop)):
                spectrum, lab = self.phrases[i].perform(label=label, 
                    keySignature=keySignature)
                return spectrum, lab
    
class Score( object ):
    """
        This Score object contains many melodies.
    """
    def __init__(self, melodies = [], directory=None, observed=None, suffix='raw'):
        self.syntheticMelodies = []
        self.syntheticMelodies = melodies
        self.directory = directory
        self.observed = observed
        self.suffix = suffix
        self.loadMelodies()

    def loadMelodies(self):
        melodyFiles = glob.glob(self.directory+'*'+self.suffix+'.fits')
        self.syntheticMelodies = []
        for melody in melodyFiles:
            print("%s" % melody)
            self.syntheticMelodies.append(SyntheticMelody(filename=melody))

        if not(self.observed==None):
            self.ObservedMelodies = [ObservedMelody.fromFile(filename=self.observed, 
                label='TWHydra')]
            self.ObservedMelodies[0].loadData()

    def getMelodyParams(self):
        raw_labels = []
        interpolated_labels = []
        integrated_labels = []
        convolved_labels = []
        observed_labels = []
        for melody in self.syntheticMelodies:
            for phrase in melody.phrases:
                for raw in phrase.rawLabels:
                    raw_labels.append(raw)
                for interpolated in phrase.interpolatedLabels:
                    interpolated_labels.append(interpolated)
                for integrated in phrase.integratedLabels:
                    integrated_labels.append(integrated)
                for convolved in phrase.convolvedLabels:
                    convolved_labels.append(convolved)
            #raw_labels.append(melody.rawLabel)
            #if len(melody.interpolatedLabels) > 0:
            #    for interpolated in melody.interpolatedLabels:
            #        interpolated_labels.append(interpolated)
            #if len(melody.integratedLabels) > 0:
            #    for integrated in melody.integratedLabels:
            #        integrated_labels.append(integrated)
            #if len(melody.convolvedLabels) > 0:
            #    for convolved in melody.convolvedLabels:
            #        convolved_labels.append(convolved)
        if not(self.observed==None):
            for melody in self.ObservedMelodies:
                for phrase in melody.phrases:
                    for obs in phrase.observedLabels:
                        observed_labels.append(obs)

        return raw_labels, interpolated_labels, integrated_labels, convolved_labels, observed_labels

    def selectEnsemble(self, selectedLabels=[], keySignature='CONVOLVED'):
        for melody in self.syntheticMelodies:
            for phrase in numpy.array(melody.phrases)[numpy.array(melody.selectedPhrases) == True]:
                if keySignature == 'RAW':
                    for raw in phrase.rawLabels:
                        if raw in selectedLabels:
                            raw.parameters["SELECTED"] = True
                        else:
                            raw.parameters["SELECTED"] = False

                if keySignature == 'INTERPOLATED':
                    for interpolated in phrase.interpolatedLabels:
                        if interpolated in selectedLabels:
                            interpolated.parameters["SELECTED"] = True
                        else:
                            interpolated.parameters["SELECTED"] = False

                if keySignature == 'INTEGRATED':
                    for integrated in phrase.integratedLabels:
                        if integrated in selectedLabels:
                            integrated.parameters["SELECTED"] = True
                        else:
                            integrated.parameters["SELECTED"] = False
                
                if keySignature == 'CONVOLVED':
                    for convolved in phrase.convolvedLabels:
                        if convolved in selectedLabels:
                            convolved.parameters["SELECTED"] = True
                        else:
                            convolved.parameters["SELECTED"] = False

    def addToEnsemble(self, selectedLabels = [], keySignature='CONVOLVED'):
        if keySignature == 'RAW':
            for melody in self.syntheticMelodies:
                if melody.rawLabel in selectedLabels:
                    melody.rawLabel.parameters["SELECTED"] = True
        elif keySignature == 'INTERPOLATED':
            for melody in self.syntheticMelodies:
                for interpolated in melody.interpolatedLabels:
                    if interpolated in selectedLabels:
                        interpolated.parameters["SELECTED"] = True
        elif keySignature == 'INTEGRATED':
            for melody in self.syntheticMelodies:
                for integrated in melody.integratedLabels:
                    if integrated in selectedLabels:
                        integrated.parameters["SELECTED"] = True
        elif keySignature == 'CONVOLVED':
            for melody in self.syntheticMelodies:
                for convolved in melody.convolvedLabels:
                    if convolved in selectedLabels:
                        convolved.parameters["SELECTED"] = True

    def selectMelodies(self, wlRange=[], selectAll=False):
        for melody in self.syntheticMelodies:
            melody.selectPhrases(wlRange=wlRange, selectAll=selectAll)

        if not(self.observed==None):
            for observed in self.ObservedMelodies:
                observed.selectPhrases(wlRange=wlRange, selectAll=selectAll)

    def tune(self, save=False):
        for melody in self.syntheticMelodies:
            print melody
            melody.tune(save=save)


    def rehearse(self, vsini=0.0, R=0.0, binToObserved=False):
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
            if binToObserved:
                melody.rehearse(vsini=vsini, R=R, 
                    observedWl = self.compositeObserved.wl)
            else:
                melody.rehearse(vsini=vsini, R=R, observedWl=None)

    def perform(self, selectedLabels=[], keySignature = "CONVOLVED"):
        spectra = []
        labels = []
        parameters = []
        for melody in self.syntheticMelodies:
            spectrum = []
            params = []
            for phrase in numpy.array(melody.phrases)[numpy.array(melody.selectedPhrases)==True]:
                if keySignature == "INTERPOLATED":
                    for interpolated in phrase.interpolatedLabels:
                        if interpolated in selectedLabels:
                            sp, p = melody.perform(interpolated,
                                keySignature=keySignature)
                            spectra.append(sp)
                            params.append(p)
                            labels.append(interpolated)

                if keySignature == "INTEGRATED":
                    for integrated in phrase.integratedLabels:
                        if integrated in selectedLabels:
                            sp, p = melody.perform(integrated,
                                    keySignature=keySignature)
                            spectra.append(sp)
                            params.append(p)
                            labels.append(integrated)

                if keySignature == "CONVOLVED":
                    for convolved in phrase.convolvedLabels:
                        if convolved.parameters["SELECTED"]:
                            if convolved in selectedLabels:
                                sp, p = melody.perform(convolved,
                                        keySignature=keySignature)
                                spectrum.append(sp)
                                params.append(p)
                                labels.append(convolved)
            if len(spectrum) > 0:
                spectra.append(spectrum)
                parameters.append(params)

        return spectra, parameters, labels


    def listen(self):   # gets the observed spectrum
        self.compositeObserved = []
        spectra = []
        labels = []
        for observed in self.ObservedMelodies:
            spectrum, label = observed.perform()
            for sp, l in zip(spectrum, label):
                spectra.append(sp)
                labels.append(l)
            #compositeSpectrum = None
            #for sp in spectrum:
            #    compositeSpectrum = SpectralTools.mergeSpectra(first=compositeSpectrum,
            #            second=sp)
                
        
            #self.compositeObserved.append(compositeSpectrum)
        #return self.compositeObserved, label
        return spectra, labels

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
