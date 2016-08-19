import astropy.io.fits as pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy
import os
import matplotlib.lines as Lines
import time

class Label( object ):
    """
    Label
    
    A Label contains information about a Spectrum, as well as a reference
    to the spectrum.
    
    Parameters are stored in a dictionary object called 'parameters'
    
    The reference to the Spectrum is stored in an object called 'reference'
    """
    def __init__(self, parameters, reference = None, Spectrum = None, Phrase = None,
            Melody = None, Score = None):
        self.parameters = parameters
        self.reference = reference
        self.Spectrum = Spectrum
        if self.Spectrum != None:
            self.Spectrum.addLabel(self)
        self.Phrase = Phrase
        self.Melody = Melody
        self.Score = Score


    def addReference(self, Spectrum=None, Phrase=None, Melody=None, Score=None,
            reference=None):
        """
        Label.addReference(reference)
        
        Adds/replaces a reference to the label.
        
        reference [SpectralTools.Spectrum] = object referenced by the Label
        """
        if Spectrum != None:
            self.Spectrum = Spectrum
            self.Spectrum.label = self
        if Phrase != None:
            self.Phrase = Phrase
        if Melody != None:
            self.Melody = Melody
        if Score != None:
            self.Score = Score
        if reference != None:
            self.reference = reference

    def copy(self):
        """
        Label.copy()
        
        provides a copy of the Label
        """
        
        parameters = {}
        for key in self.parameters.keys():
            parameters[key] = self.parameters[key]
        return Label(parameters, reference = self.reference, Spectrum = self.Spectrum,
                Phrase = self.Phrase, Melody = self.Melody, Score = self.Score)

    def merge(self, other):
        """
        merged = Label.merge(other)
        
        Merges the current label with another label.  The merge function
        is used when spectra from the same underlying model (Teff, log g,
        B-field, Resolving power, etc...) are merged together.  If the
        spectra referred to by the two labels do not have the same underlying
        model, a Moog960 Error is thrown.
        
        other [Label] = other label.
        """
        parameters = {}
        for key in self.parameters.keys():
            if not(key in ['WLSTART', 'WLSTOP', 'SELECTED', 'LABEL']):
                if self.parameters[key] != other.parameters[key]:
                    raise Moog960Error(2, '[%s] %s != %s' % (key, self.parameters[key],
                            other.parameters[key]))
            parameters[key] = self.parameters[key]
        oldStart = self.parameters["WLSTART"]
        oldStop = self.parameters["WLSTOP"]
        newStart = other.parameters["WLSTART"]
        newStop = other.parameters["WLSTOP"]
        parameters["WLSTART"] = numpy.min([oldStart, newStart])
        parameters["WLSTOP"] = numpy.max([oldStop, newStop])

        return Label(parameters)

    def __eq__(self, other):
        """
        Label.__eq__(other)
        
        Determines whether or not two labels are equal to one another.
        """
        if other == None:
            return False
        for key in self.parameters.keys():
            if self.parameters[key] != other.parameters[key]:
                return False
        return True

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
        self.message[2] = "Blending Error!  Requested parameter is outside grid!! %s" %errmsg
        self.message[3] = "Merging Error! The two spectra do not have the same parameters!!\n %s" % errmsg

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

    def inWlRange(self, wlStart, wlStop, exact = False):
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
        if exact:
            return ((self.wlStart == wlStart) & (self.wlStop == wlStop))
        else:
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
        parameters["SELECTED"] = False
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
    #def __init__(self, rawData=[], interpolatedData=[], 
    #             integratedData=[], convolvedData=[], 
    #             rawLabels=[], interpolatedLabels=[],
    #             integratedLabels=[], convolvedLabels=[], diskInt = None):
    def __init__(self, rawData=[], interpolatedData=[], 
                 integratedData=[], convolvedData=[], 
                 diskInt = None, parent=None):
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
        #self.rawLabels = rawLabels
        self.interpolatedData = interpolatedData
        #self.interpolatedLabels = interpolatedLabels
        self.integratedData = integratedData
        #self.integratedLabels = integratedLabels
        self.convolvedData = convolvedData
        #self.convolvedLabels = convolvedLabels
        self.parent = parent
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
                 sourceType="RAW", parameters=None, parent=None):
        """
        Creates a phrase from a data file    
        """
        parameters["PHRASE"] = self
        if sourceType =="RAW":
            rawData = []
            rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
                filename=filename, ext=ext))
            interpolatedData = []
            integratedData = []
            convolvedData = []
            
            if parent != None:
                parent.parent.raw_labels.append(Label(parameters, reference=rawData[-1], Phrase=self))
            #rawLabels = []
            #rawLabels.append(Label(parameters, reference=rawData[-1], Phrase=self))
            #interpolatedLabels = []
            #integratedLabels = []
            #convolvedLabels = []
        elif sourceType == "INTERPOLATED":
            rawData = []
            interpolatedData = []
            interpolatedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            integratedData = []
            convolvedData = []

            if parent != None:
                parent.parent.interpolated_labels.append(Label(parameters, reference=interpolatedData[-1], Phrase=self))
            #rawLabels = []
            #interpolatedLabels = []
            #interpolatedLabels.append(Label(parameters, reference = interpolatedData[-1],
            #    Phrase=self))
            #integratedLabels = []
            #convolvedLabels = []
        elif sourceType == "INTEGRATED":
            rawData = []
            interpolatedData = []
            integratedData = []
            integratedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            convolvedData = []

            if parent != None:
                parent.parent.integrated_labels.append(Label(parameters, reference=integratedData[-1], Phrase=self))
            
            #rawLabels = []
            #interpolatedLabels = []
            #integratedLabels = []
            #integratedLabels.append(Label(parameters, reference = integratedData[-1],
            #    Phrase=self))
            #convolvedLabels = []
        elif sourceType == "CONVOLVED":
            rawData = []
            interpolatedData = []
            integratedData = []
            convolvedData = []
            convolvedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
                

            if parent != None:
                parent.parent.convolved_labels.append(Label(parameters, reference=convolvedData[-1], Phrase=self,
                       Spectrum=convolvedData[-1]))

            #rawLabels = []
            #interpolatedLabels = []
            #integratedLabels = []
            #convolvedLabels = []
            #convolvedLabels.append(Label(parameters, reference=convolvedData[-1], 
            #    Phrase=self))

        #return self(rawData=rawData, interpolatedData=interpolatedData, 
        #        integratedData=integratedData, convolvedData=convolvedData, 
        #        rawLabels=rawLabels, interpolatedLabels=interpolatedLabels,
        #        integratedLabels=integratedLabels, convolvedLabels=convolvedLabels, diskInt=diskInt)
        return self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                diskInt=diskInt, parent=parent)
    
    def addSpectrum(self, hdr, data=None, filename=None, ext=None, sourceType="RAW", parameters=None):
        """
        Moog960::SyntheticPhrase::addSpectrum(hdr, data=None, filename=None, ext=None, sourceType="RAW")
        
        Adds a spectrum to the phrase, depending on its source (RAW, INTERPOLATED, INTEGRATED, CONVOLVED)
        """
        if sourceType=="RAW":
            self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.rawLabels.append(Label(parameters, reference=self.rawData[-1]))
        elif sourceType =="INTERPOLATED":
            self.processedData.interpolated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.interpolatedLabels.append(Label(parameters, 
                reference=self.processedData.interpolated[-1]))
        elif sourceType =="INTEGRATED":
            self.processedData.integrated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.integratedLabels.append(Label(parameters, 
                reference = self.processedData.integrated[-1]))
        elif sourceType =="CONVOLVED":
            self.processedData.convolved.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.convolvedLabels.append(Label(parameters, 
                reference = self.processedData.convolved[-1]))

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
            self.integratedLabels.append(Label(parameters, 
                reference = self.processedData.integrated[-1]))
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
            parameters["SELECTED"] = False
            parameters["LABEL"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini)
            self.integratedLabels.append(Label(parameters,
                reference = self.processedData.integrated[-1]))
        if 'CONVOLVED' in created:
            parameters = self.rawLabels[0].parameters.copy()
            parameters["VSINI"] = vsini
            parameters["R"] = R
            parameters["WLSTART"] = self.wlStart
            parameters["WLSTOP"] = self.wlStop
            parameters["SELECTED"] = False
            parameters["LABEL"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini, R)

            self.convolvedLabels.append(Label(parameters,
                reference = self.processedData.convolved[-1]))

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
        
        """
        if keySignature == "CONVOLVED":
            if label == None:
                ret_labels = []
                for convolved in self.convolvedData:
                    print "Blah"
                    raw_input()
                    ret_labels.append(convolved.label)
                return self.convolvedData, ret_labels
            else:
                for l, convolved in zip(self.parent.parent.convolved_labels, self.convolvedData):
                    if l==label:
                        return convolved, label
        #"""
        
        return label.Spectrum, label
        
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
            SpectrumHDU.name = "%.4fA - %.4fA PHI=%.3f MU=%.3f DELTAV=%.3f" % \
                     (hdr.get("wlStart"), hdr.get("wlStop"), hdr.get("PHI_ANGLE"), hdr.get("MU"), hdr.get('DELTAV'))

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
                #index = self.convolvedLabels.index(label)
                spectrum = label.Spectrum
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
    def __init__(self, phrases = [], filename=None, label=None, header=None, Score=None):
        self.phrases = phrases
        self.header = header
        self.nPhrases = len(self.phrases)
        self.filename = filename
        self.selectedPhrases = [False for i in range(self.nPhrases)]
        self.muted = True
        self.label = label
        self.Score = Score   # Pointer to Score which contains this melody (optional)
        
    def addPhrases(self, phrases=[]):    # WTF?  This doesn't look right
        for phrase in phrases:
            self.phrases.append(phrase)
        self.nPhrases = len(self.phrases)

    def selectPhrases(self, wlRange=[], exact=False, selectAll=False):
        """
        Moog960::Melody::selectPhrases(wlRange=[], selectAll=False):
        
        selectPhrases() goes through the phrases in a melody and selects
            those phrases which overlap with the set wavelength range
            
        Input:
            wlRange = [wlStart, wlStart]
            
        Output:
            Nothing
        """

        for phrase in self.phrases:
            if selectAll:
                for convolved in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = True
            else:
                for convolved in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = phrase.inWlRange(wlStart=wlRange[0],
                        wlStop=wlRange[1], exact=exact)

    def record(self, labels=None, basename = None):
        """
        Moog960::Melody::record(self, labels=[], basename=None)
        
        record() saves the selected phrases of a melody to disk
        """
        
        if self.filename != None:
            filename = self.filename[:-4]+'_saved.fits'
        else:
            filename = 'RecordedSpectrum.fits'
        if labels != None:
            for label in labels:
                R = label.parameters["R"]
                vsini = label.parameters["VSINI"]
                if basename != None:
                    filename = basename+"_T%d_G%.2f_B%.2f_R%d_V%.2f.fits" % (
                        self.Teff, self.logg, self.B, R, vsini)
                print("Recording record \'%s\' to disk" % label.parameters["LABEL"])
                label.parameters["PHRASE"].saveConvolved(label=label)
                #for i in range(self.nPhrases):
                #    if self.selectedPhrases[i]:
                #       self.phrases[i].saveConvolved(vsini=vsini, R=R,
                #               filename=filename, header=self.header, 
                #               wlStart=label.parameters["WLSTART"], 
                #               wlStop = label.parameters["WLSTOP"])
        else:
            for i in range(self.nPhrases):
                if self.selectedPhrases[i]:
                    self.phrases[i].saveConvolved()

                           
            #TODO make Labels actually do what they are supposed to.


class ObservedMelody( Melody ):
    def __init__(self, phrases = [], filename=None, label=None, header=None, Score=None):
        super(ObservedMelody, self).__init__(phrases=phrases, filename=filename, 
                label=label, header=header, Score=Score)

    @classmethod
    def fromFile(self, filename=None, label=None, parent=None):
        info = pyfits.info(filename, output='')
        nPhrases = len(info)-1
        header = pyfits.getheader(filename, ext=0)
        phrases = []
        for i in range(nPhrases):
            hdr = pyfits.getheader(filename, ext=i+1)
            phrases.append(ObservedPhrase.fromFile(header=hdr, filename=filename, ext=i+1))

        return self(phrases=phrases, filename=filename, label=label, header=header, parent=parent)
        
    def loadData(self):
        for phrase in self.phrases:
            phrase.loadData()
            
    def rehearse(self, **kwargs):
        return

    def selectPhrases(self, wlRange=[], exact=False, selectAll=False):
        self.selectedPhrases = []
        for phrase in self.phrases:
            if selectAll:
                self.selectedPhrases.append(True)
            else:
                self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
                    wlStop=wlRange[1], exact=exact))

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
    def __init__(self, phrases = [], filename=None, label=None, header=None, Score=None):
        super(SyntheticMelody, self).__init__(phrases = phrases, filename=filename,
                label=label, header=header, Score=Score)
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
        parameters["SELECTED"] = False

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
                    sourceType=self.contents, parameters=parameters, parent=self))


        self.nPhrases = len(self.phrases)

    def tune(self, save=False):
        for i in range(self.nPhrases):
            self.phrases[i].tune(save=save, header=self.header)

    def selectPhrases(self, wlRange=[], exact=False, selectAll=False, keySignature="CONVOLVED"):
        for phrase in self.phrases:
            if selectAll:
                for convolved in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = True
            else:
                for convolved in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = phrase.inWlRange(wlStart=wlRange[0],
                        wlStop=wlRange[1], exact=exact)
        #self.selectedPhrases = []
        #for phrase in self.phrases:
        #    if selectAll:
        #        self.selectedPhrases.append(True)
        #    else:
        #        self.selectedPhrases.append(phrase.inWlRange(wlStart=wlRange[0],
        #            wlStop=wlRange[1], exact=exact))

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
            if (self.phrases[i].convolvedData[0].label.parameters["SELECTED"] &
                    (label.parameters["WLSTART"] == self.phrases[i].wlStart) &
                    (label.parameters["WLSTOP"] == self.phrases[i].wlStop)):
                spectrum, lab = self.phrases[i].perform(label=label, 
                    keySignature=keySignature)
                return spectrum, lab

    def master(self, keySignature='CONVOLVED'):
        spectra = []
        labels = []
        for i in range(self.nPhrases):
             spectra.append(self.phrases[i].convolvedData[0])
             labels.append(self.phrases[i].convolvedData[0].label)
        spectra = numpy.array(spectra).T
        labels = numpy.array(labels).T

        order = numpy.argsort(labels)
        spectra = spectra[order]
        labels = labels[order]

        mergedSpectra = spectra[0]
        mergedLabel = labels[0]
        for sp, l in zip(spectra[1:], labels[1:]):
            mergedSpectra = mergedSpectra.mergeSpectra(sp)
            mergedLabel = mergedLabel.merge(l)

        mergedLabel.parameters["SELECTED"] = True
        mergedLabel.addReference(Spectrum = mergedSpectra)
        if keySignature=='CONVOLVED':
            newPhrase = SyntheticPhrase(convolvedData=[mergedSpectra], parent=self)
            self.parent.convolved_labels.append(mergedLabel)
            self.addPhrases(phrases = [newPhrase])
        
        return mergedLabel
        

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
        self.computeGridPoints()
        #self.getMelodyParams(retLabels = False)

    def loadMelodies(self):
        melodyFiles = glob.glob(self.directory+'*'+self.suffix+'.fits')
        self.syntheticMelodies = []
        self.raw_labels = []
        self.interpolated_labels = []
        self.integrated_labels = []
        self.convolved_labels = []
        self.observed_labels = []
        for melody in melodyFiles:
            print("%s" % melody)
            self.syntheticMelodies.append(SyntheticMelody(filename=melody, parent=self))

        if not(self.observed==None):
            self.ObservedMelodies = [ObservedMelody.fromFile(filename=self.observed, 
                label='TWHydra', parent=self)]
            self.ObservedMelodies[0].loadData()

    def getMelodyParams(self, retLabels=True):
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
        if not(self.observed==None):
            for melody in self.ObservedMelodies:
                for phrase in melody.phrases:
                    for obs in phrase.observedLabels:
                        observed_labels.append(obs)

        if retLabels:
            return raw_labels, interpolated_labels, integrated_labels, convolved_labels, observed_labels
        else:
            self.raw_labels = raw_labels
            self.interpolated_labels = interpolated_labels
            self.integrated_labels = integrated_labels
            self.convolved_labels = convolved_labels
            self.observed_labels = observed_labels
            
            # raw GridPoints
            params = {}
            self.rawGridPoints = {}
            keys = ["TEFF", "LOGG", "BFIELD"]
            for key in keys:
                params[key] = []
            for raw in self.raw_labels:
                for key in raw.keys():
                    params[key].append(raw[key])
            for key in params.keys():
                self.rawGridPoints[key] = numpy.unique(params[key])

            # interpolated GridPoints
            params = {}
            self.interpolatedGridPoints = {}
            keys = ["TEFF", "LOGG", "BFIELD"]
            for key in keys:
                params[key] = []
            for interpolated in self.interpolated_labels:
                for key in interpolated.keys():
                    params[key].append(interpolated[key])
            for key in params.keys():
                self.interpolatedGridPoints[key] = numpy.unique(params[key])

            # integrated GridPoints
            params = {}
            self.integratedGridPoints = {}
            keys = ["TEFF", "LOGG", "BFIELD"]
            for key in keys:
                params[key] = []
            for integrated in self.integrated_labels:
                for key in integrated.keys():
                    params[key].append(integrated[key])
            for key in params.keys():
                self.integratedGridPoints[key] = numpy.unique(params[key])

            # convolved GridPoints
            params = {}
            self.convolvedGridPoints = {}
            keys = ["TEFF", "LOGG", "BFIELD", "WLSTART", "WLSTOP", "SELECTED", "VSINI", 
                    "LABEL", "R"]
            for key in keys:
                params[key] = []
            for convolved in self.convolved_labels:
                for key in keys:#convolved.parameters.keys():
                    params[key].append(convolved.parameters[key])
            for key in params.keys():
                self.convolvedGridPoints[key] = numpy.unique(params[key])

            # observed GridPoints
            params = {}
            self.observedGridPoints = {}
            keys = ["WLSTART", "WLSTOP"]
            for key in keys:
                params[key] = []
            for observed in self.observed_labels:
                for key in observed.parameters.keys():
                    params[key].append(observed.parameters[key])
            for key in params.keys():
                self.observedGridPoints[key] = numpy.unique(params[key])

    def computeGridPoints(self):
        # raw GridPoints
        params = {}
        self.rawGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD"]
        for key in keys:
            params[key] = []
        for raw in self.raw_labels:
            for key in raw.keys():
                params[key].append(raw[key])
        for key in params.keys():
            self.rawGridPoints[key] = numpy.unique(params[key])

        # interpolated GridPoints
        params = {}
        self.interpolatedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD"]
        for key in keys:
            params[key] = []
        for interpolated in self.interpolated_labels:
            for key in interpolated.keys():
                params[key].append(interpolated[key])
        for key in params.keys():
            self.interpolatedGridPoints[key] = numpy.unique(params[key])

        # integrated GridPoints
        params = {}
        self.integratedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD"]
        for key in keys:
            params[key] = []
        for integrated in self.integrated_labels:
            for key in integrated.keys():
                params[key].append(integrated[key])
        for key in params.keys():
            self.integratedGridPoints[key] = numpy.unique(params[key])

        # convolved GridPoints
        params = {}
        self.convolvedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD", "WLSTART", "WLSTOP", "SELECTED", "VSINI", 
                "LABEL", "R"]
        for key in keys:
            params[key] = []
        for convolved in self.convolved_labels:
            for key in keys:#convolved.parameters.keys():
                params[key].append(convolved.parameters[key])
        for key in params.keys():
            self.convolvedGridPoints[key] = numpy.unique(params[key])

        # observed GridPoints
        params = {}
        self.observedGridPoints = {}
        keys = ["WLSTART", "WLSTOP"]
        for key in keys:
            params[key] = []
        for observed in self.observed_labels:
            for key in observed.parameters.keys():
                params[key].append(observed.parameters[key])
        for key in params.keys():
            self.observedGridPoints[key] = numpy.unique(params[key])

    def selectEnsemble(self, selectedLabels=[], keySignature='CONVOLVED'):
        for melody in self.syntheticMelodies:
            for phrase in numpy.array(melody.phrases)[numpy.array(melody.selectedPhrases) == True]:
                if keySignature == 'RAW':
                    for raw in self.raw_labels:
                        if raw in selectedLabels:
                            raw.parameters["SELECTED"] = True
                        else:
                            raw.parameters["SELECTED"] = False

                if keySignature == 'INTERPOLATED':
                    for interpolated in self.interpolated_labels:
                        if interpolated in selectedLabels:
                            interpolated.parameters["SELECTED"] = True
                        else:
                            interpolated.parameters["SELECTED"] = False

                if keySignature == 'INTEGRATED':
                    for integrated in self.integrated_labels:
                        if integrated in selectedLabels:
                            integrated.parameters["SELECTED"] = True
                        else:
                            integrated.parameters["SELECTED"] = False
                
                if keySignature == 'CONVOLVED':
                    for convolved in self.convolved_labels:
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

    def selectExcerpt(self, wlRange=[], exact=False, selectAll=False):
        for melody in self.syntheticMelodies:
            melody.selectPhrases(wlRange=wlRange, exact=exact, selectAll=selectAll)

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
            if keySignature == "CONVOLVED":
                for convolved in self.convolved_labels:
                    if convolved.parameters["SELECTED"]:
                        if convolved in selectedLabels:
                            
                            #sp, p = melody.perform(convolved,
                            #        keySignature=keySignature)
                            spectrum.append(convolved.Spectrum)
                            params.append(convolved)
                            labels.append(convolved)
            """
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
                    for convolved in self.parent.convolved_labels:
                        if convolved.parameters["SELECTED"]:
                            if convolved in selectedLabels:
                                sp, p = melody.perform(convolved,
                                        keySignature=keySignature)
                                spectrum.append(sp)
                                params.append(p)
                                labels.append(convolved)

                if keySignature == "MERGED":
                    for merged in phrase.mergedLabels:
                        if merged.parameters["SELECTED"]:
                            if merged in mergedLabels:
                                sp, p = melody.perform(merged,
                                        keySignature=keySignature)
                                spectrum.append(sp)
                                params.append(p)
                                labels.append(merged)
            #"""

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

    def master(self, selectedLabels=[], keySignature='CONVOLVED'):
        """
           Moog960::SyntheticMelody::master()
           This routine merges all phrases into a single master phrase
        """
        #compositePhrase = SpectralTools.Spectrum()
        #self.selectMelodies(selectAll=True)
        #self.selectEnsemble(selectedLabels=selectedLabels, keySignature=keySignature)
        mastered = []
        for synthetic in self.syntheticMelodies:
            mastered.append(synthetic.master(keySignature=keySignature))

        return mastered

        #self.getMelodyParams(retLabels=False)
        #spectra, parameters, labels = self.perform(selectedLabels=selectedLabels, 
        #        keySignature='CONVOLVED')


    def blend(self, desiredParameters={}):
        gridPoints = {}
        for key in desiredParameters.keys():
            points = self.convolvedGridPoints[key]
            low = numpy.min(points)
            high = numpy.max(points)
            desired = desiredParameters[key]
            print desired
            if ( (desired > high) | (desired < low) ):
                raise Moog960Error (2, "%.2f not within range (%.2f, %.2f)" % (desired, low, high))
            gridPoints[key] = [numpy.sort(points[points<=desired])[-1],
                                    numpy.sort(points[points >=desired])[0]]

        selected = []
        for label in self.convolved_labels:
            if label.parameters["TEFF"] in gridPoints["TEFF"]:
                if label.parameters["LOGG"] in gridPoints["LOGG"]:
                    if label.parameters["BFIELD"] in gridPoints["BFIELD"]:
                        if label.parameters["SELECTED"]:
                            selected.append(label)

        spectra = []
        params = []
        for label in selected:
            spectra.append(label.Spectrum)
            params.append(label)
            

        Bspectra = []
        Bparams = []
        if len(params) == 1:
            Bspectra.append(spectra[0])
            Bparams.append(params[0])
        else:
            while len(params) > 0:
                sp1 = spectra.pop(0)
                p1 = params.pop(0)
                p2 = p1
                for i in range(len(params)):
                    if ( (params[i].parameters["TEFF"] == p1.parameters["TEFF"]) &
                            (params[i].parameters["LOGG"] == p1.parameters["LOGG"])):
                        sp2 = spectra.pop(i)
                        p2 = params.pop(i)
                        break

                if p2 != p1:
                    newp = []
                    fraction = -1.0
                    order1 = numpy.argsort(numpy.array(p1))
                    order2 = numpy.argsort(numpy.array(p2))
                    p1 = numpy.array(p1)[order1].tolist()
                    p2 = numpy.array(p2)[order2].tolist()
                    sp1 = numpy.array(sp1)[order1].tolist()
                    sp2 = numpy.array(sp2)[order2].tolist()
                    for param1, param2 in zip(p1, p2):
                        newp = p1.copy()
                        for key in desiredParameters.keys():
                            if param1.parameters[key] != param2.parameters[key]:
                                distance = param1.parameters[key] - param2.parameters[key]
                                fraction = (param1.parameters[key] - desiredParameters[key])/distance
                                newp.parameters[key] = desiredParameters[key]
                        Bparams.append(newp)
                    if fraction != -1.0:
                        for s1, s2 in zip(sp1, sp2):
                            Bspectra.append(s1.blend(s2, fraction))
                    else:
                        for s1 in sp1:
                            Bspectra.append(s1)

                else:
                    Bspectra.append(sp1)
                    Bparams.append(p1)

        Gspectra = []
        Gparams = []
        if len(Bparams) == 1:
            Gspectra.append(Bspectra[0])
            Gparams.append(Bparams[0])
        else:
            while len(Bparams) > 0:
                sp1 = Bspectra.pop(0)
                p1 = Bparams.pop(0)
                p2 = p1
                for i in range(len(Bparams)):
                    if ( (Bparams[i].parameters["TEFF"] == p1.parameters["TEFF"])):
                        sp2 = Bspectra.pop(i)
                        p2 = []
                        p2 = Bparams.pop(i)
                        break
                if p2 != p1:
                    newp = []
                    fraction = -1.0
                    order1 = numpy.argsort(numpy.array(p1))
                    order2 = numpy.argsort(numpy.array(p2))
                    p1 = numpy.array(p1)[order1].tolist()
                    p2 = numpy.array(p2)[order2].tolist()
                    sp1 = numpy.array(sp1)[order1].tolist()
                    sp2 = numpy.array(sp2)[order2].tolist()
                    for param1, param2 in zip(p1, p2):
                        newp = p1.copy()
                        for key in desiredParameters.keys():
                            if param1.parameters[key] != param2.parameters[key]:
                                distance = param1.parameters[key] - param2.parameters[key]
                                fraction = (param1.parameters[key] - desiredParameters[key])/distance
                                newp.parameters[key] = desiredParameters[key]
                        Gparams.append(newp)
                    if fraction != -1.0:
                        for s1, s2 in zip(sp1, sp2):
                            Gspectra.append(s1.blend(s2, fraction))
                    else:
                        for s1 in sp1:
                            Gspectra.append(s1)
                else:
                    Gspectra.append(sp1)
                    Gparams.append(p1)

        interpolatedSpectrum = []
        interpolatedParams = []

        if len(Gparams) == 1:
            interpolatedSpectrum = [Gspectra[0]]
            interpolatedParams = [Gparams[0]]
        else:
            newp = []
            fraction = -1.0
            distance = Gparams[0].parameters["TEFF"] - Gparams[1].parameters["TEFF"]
            fraction = (Gparams[0].parameters["TEFF"] - desiredParameters["TEFF"])/distance
            newp = Gparams[0].copy()
            newp.parameters["TEFF"] = desiredParameters["TEFF"]
            interpolatedParams.append(newp)
            #order1 = numpy.argsort(numpy.array(Gparams[0]))
            #order2 = numpy.argsort(numpy.array(Gparams[1]))
            #p1 = numpy.array(Gparams[0])[order1].tolist()
            #p2 = numpy.array(Gparams[1])[order2].tolist()
            #sp1 = numpy.array(Gspectra[0])[order1].tolist()
            #sp2 = numpy.array(Gspectra[1])[order2].tolist()
            #for param1, param2 in zip(p1, p2):
            #    newp = param1.copy()
            #    for key in desiredParameters.keys():
            #        if param1.parameters[key] != param2.parameters[key]:
            #            distance = param1.parameters[key] - param2.parameters[key]
            #            fraction = (param1.parameters[key] - desiredParameters[key])/distance
            #            newp.parameters[key] = desiredParameters[key]
            #    interpolatedParams.append(newp)
            if fraction != -1.0:
                #for s1, s2 in zip(sp1, sp2):
                #    interpolatedSpectrum.append(s1.blend(s2, fraction))
                interpolatedSpectrum.append(Gspectra[0].blend(Gspectra[1], fraction))
            else:
                for s1 in sp1:
                    interpolatedSpectrum.append(s1)

        phrases = []
        for sp, p in zip(interpolatedSpectrum, interpolatedParams):
            p.addReference(sp)
            self.convolved_labels.append(p)
            phrases.append(SyntheticPhrase(convolvedData=[sp], diskInt='BEACHBALL'))
        header=pyfits.Header()
        for key in desiredParameters.keys():
            header.set(key, desiredParameters[key])
        header.set("SPECTRUM_CONTENTS", "CONVOLVED")
        blended = SyntheticMelody(phrases=phrases, header=header, parent=self)
        #blended.selectPhrases(selectAll=True)
        #print len(interpolatedParams)
        #raw_input()
        #self.selectMelodies(selectAll=True)

        return blended, interpolatedParams

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
