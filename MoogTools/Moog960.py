import astropy.io.fits as pyfits
import AstroUtils
import glob
import MoogTools
import SpectralTools
import numpy
import os
import matplotlib.lines as Lines
import time
import random
import string

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
            if not(key in ['WLSTART', 'WLSTOP', 'SELECTED', 'NAME', 'PHRASE']):
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

    def __le__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() <= (other.getWlStart())

    def __lt__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() < (other.getWlStart())

    def __eq__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() == (other.getWlStart())

    def __ne__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() != (other.getWlStart())

    def __ge__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() >= (other.getWlStart())

    def __gt__(self, other):
        if hasattr(other, 'getWlStart'):
            return self.getWlStart() > (other.getWlStart())

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
    def __init__(self, wlStart, wlStop, ID=None):
        self.wlStart = wlStart
        self.wlStop = wlStop
        if ID == None:
            self.ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        else:
           self.ID = ID
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
    def __init__(self, observedData = [], Melody=None, ID=None):
        wlStart = observedData[0].header.get('WLSTART')
        wlStop = observedData[0].header.get('WLSTOP')
        super(ObservedPhrase, self).__init__(wlStart, wlStop, ID=ID)

        self.observedData = observedData
        #self.observedLabels = observedLabels

    @classmethod
    def fromFile(self, header=None, data=None, filename=None, ext=None, parameters=None, Melody=None):
        observed = SpectralTools.Spectrum.from_file(header= header, data=data,
                filename=filename, ext=ext)
        parameters = {}
        ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        parameters["PHRASE"] = ID
        parameters["MELODY"] = Melody.ID
        parameters["SCORE"] = Melody.Score.ID
        parameters["WLSTART"] = observed.header.get('WLSTART')
        parameters["WLSTOP"] = observed.header.get('WLSTOP')
        parameters["SELECTED"] = False

        retval =  self(observedData=[observed], Melody=Melody, ID=ID)
        if Melody != None:
            if not(Melody.ID in Melody.Score.observed_labels.keys()):
                Melody.Score.observed_labels[Melody.ID] = {}
            Melody.Score.observed_labels[Melody.ID][ID] = []
            Melody.Score.observed_labels[Melody.ID][ID].append(Label(parameters, reference=retval.observedData[-1], Phrase=retval, Spectrum = retval.observedData[-1]))
        return retval
        
    def loadData(self):
        for observed in self.observedData:
            observed.loadData()
        
    def listen(self):
        return self.observedData
        
    def record(self, filename = None, dI=False, primaryHeaderKWs=None):
        self.save(filename=filename, dI=dI, primaryHeaderKWs=primaryHeaderKWs)

    def save(self, filename = None, dI=False, primaryHeaderKWs=None):
        HDUs = []
        for spectrum in self.observedData:
            hdr = spectrum.header.copy()
            spectrum.preserve(continuum=False, V=False, dI=dI)
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
            if primaryHeaderKWs != None:
                for key in primaryHeaderKWs.keys():
                    primary.header.set(key, primaryHeaderKWs[key])
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
                 diskInt = None, Melody=None, ID=None):
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
        self.interpolatedData = interpolatedData
        self.integratedData = integratedData
        self.convolvedData = convolvedData
        self.Melody = Melody
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
        if ID == None:
            self.ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        else:
           self.ID = ID

    @classmethod
    def fromFile(self, hdr, data=None, filename=None, ext=None, diskInt=None,
                 sourceType="RAW", parameters=None, Melody=None):
        """
        Creates a phrase from a data file    
        """
        ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        parameters["PHRASE"] = ID
        parameters["MELODY"] = Melody.ID
        try:
            parameters["SCORE"] = Melody.Score.ID
        except:
            parameters["SCORE"] = 1
        if sourceType =="RAW":
            rawData = []
            rawData.append(SpectralTools.Spectrum.from_file(hdr,data=data,
                filename=filename, ext=ext))
            interpolatedData = []
            integratedData = []
            convolvedData = []
            retval =  self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                diskInt=diskInt, Melody=Melody, ID=ID)
            
            if Melody != None:
                if not(Melody.ID in Melody.Score.raw_labels.keys()):
                    Melody.Score.raw_labels[Melody.ID] = {}
                Melody.Score.raw_labels[Melody.ID][ID] = []
                Melody.Score.raw_labels[Melody.ID][ID].append(Label(parameters, reference=rawData[-1], Phrase=retval))
        elif sourceType == "INTERPOLATED":
            rawData = []
            interpolatedData = []
            interpolatedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            integratedData = []
            convolvedData = []
            retval =  self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                diskInt=diskInt, Melody=Melody, ID=ID)

            if Melody != None:
                if not(Melody.ID in Melody.Score.interpolated_labels.keys()):
                    Melody.Score.interpolated_labels[Melody.ID] = {}
                Melody.Score.interpolated_labels[Melody.ID][ID] = []
                Melody.Score.interpolated_labels[Melody.ID][ID].append(Label(parameters, reference=interpolatedData[-1], Phrase=retval))
        elif sourceType == "INTEGRATED":
            rawData = []
            interpolatedData = []
            integratedData = []
            integratedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            convolvedData = []
            retval =  self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                diskInt=diskInt, Melody=Melody, ID=ID)

            if Melody != None:
                if not(Melody.ID in Melody.Score.integrated_labels.keys()):
                    Melody.Score.integrated_labels[Melody.ID] = {}
                Melody.Score.integrated_labels[Melody.ID][ID] = []
                Melody.Score.integrated_labels[Melody.ID][ID].append(Label(parameters, reference=integratedData[-1], Phrase=retval))
            
        elif sourceType == "CONVOLVED":
            rawData = []
            interpolatedData = []
            integratedData = []
            convolvedData = []
            convolvedData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            retval =  self(rawData=rawData, interpolatedData=interpolatedData, 
                integratedData=integratedData, convolvedData=convolvedData, 
                diskInt=diskInt, Melody=Melody, ID=ID)
                

            if Melody != None:
                if not(Melody.ID in Melody.Score.convolved_labels.keys()):
                    Melody.Score.convolved_labels[Melody.ID] = {}
                Melody.Score.convolved_labels[Melody.ID][ID] = []
                Melody.Score.convolved_labels[Melody.ID][ID].append(Label(parameters, reference=convolvedData[-1], Phrase=retval,
                       Spectrum=convolvedData[-1]))

        return retval
    
    def addSpectrum(self, hdr, data=None, filename=None, ext=None, sourceType="RAW", parameters=None):
        """
        Moog960::SyntheticPhrase::addSpectrum(hdr, data=None, filename=None, ext=None, sourceType="RAW")
        
        Adds a spectrum to the phrase, depending on its source (RAW, INTERPOLATED, INTEGRATED, CONVOLVED)
        """
        if sourceType=="RAW":
            self.rawData.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.Melody.Score.raw_labels[self.Melody.ID][self.ID].append(Label(parameters,
                reference=self.rawData[-1], 
                Score=self.Melody.Score, Melody=self.Melody, Spectrum=self.rawData[-1], Phrase=self))
        elif sourceType =="INTERPOLATED":
            self.processedData.interpolated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.Melody.Score.interpolated_labels[self.Melody.ID][self.ID].append(Label(parameters, 
                reference=self.processedData.interpolated[-1], Score=self.Melody.Score,
                Melody=self.Melody, Spectrum=self.processedData.interpolated[-1], Phrase=self))
        elif sourceType =="INTEGRATED":
            self.processedData.integrated.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.Melody.Score.integrated_labels[self.Melody.ID][self.ID].append(Label(parameters, 
                reference = self.processedData.integrated[-1], Score=self.Melody.Score,
                Melody=self.Melody, Spectrum=self.processedData.integrated[-1], Phrase=self))
        elif sourceType =="CONVOLVED":
            self.processedData.convolved.append(SpectralTools.Spectrum.from_file(hdr, data=data,
                filename=filename, ext=ext))
            self.Melody.Score.convolved_labels[self.Melody.ID][self.ID].append(Label(parameters, 
                reference = self.processedData.convolved[-1], Score=self.Melody.Score,
                Melody=self.Melody, Spectrum=self.processedData.convolved[-1], Phrase=self))

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
            parameters = self.Melody.Score.raw_labels[self.Melody.ID][self.ID][0].parameters.copy()
            parameters["VSINI"] = vsini
            parameters["WLSTART"] = self.wlStart
            parameters["WLSTOP"] = self.wlStop
            parameters["SELECTED"] = False
            parameters["NAME"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini)
            if self.Melody.ID in self.Melody.Score.integrated_labels.keys():
                if not(self.ID in self.Melody.Score.integrated_labels[self.Melody.ID].keys()):
                    self.Melody.Score.integrated_labels[self.Melody.ID][self.ID] = []
                self.Melody.Score.integrated_labels[self.Melody.ID][self.ID].append(Label(parameters,
                     reference = self.processedData.integrated[-1], Phrase=self, 
                     Spectrum=self.processedData.integrated[-1], Melody=self.Melody, 
                     Score=self.Melody.Score))
            else:
                self.Melody.Score.integrated_labels[self.Melody.ID] = {}
                self.Melody.Score.integrated_labels[self.Melody.ID][self.ID] = []
                self.Melody.Score.integrated_labels[self.Melody.ID][self.ID].append(Label(parameters,
                     reference = self.processedData.integrated[-1], Phrase=self, 
                     Spectrum=self.processedData.integrated[-1], Melody=self.Melody, 
                     Score=self.Melody.Score))

        if 'CONVOLVED' in created:
            parameters = self.Melody.Score.raw_labels[self.Melody.ID][self.ID][0].parameters.copy()
            parameters["VSINI"] = vsini
            parameters["R"] = R
            parameters["WLSTART"] = self.wlStart
            parameters["WLSTOP"] = self.wlStop
            parameters["SELECTED"] = False
            parameters["NAME"] = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % ( 
                            parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"],
                            vsini, R)

            if self.Melody.ID in self.Melody.Score.convolved_labels.keys():
                if not(self.ID in self.Melody.Score.convolved_labels[self.Melody.ID].keys()):
                    self.Melody.Score.convolved_labels[self.Melody.ID][self.ID] = []
                self.Melody.Score.convolved_labels[self.Melody.ID][self.ID].append(Label(parameters,
                     reference = self.processedData.convolved[-1], Phrase=self, 
                     Spectrum=self.processedData.convolved[-1], Melody=self.Melody, 
                     Score=self.Melody.Score))
            else:
                self.Melody.Score.convolved_labels[self.Melody.ID] = {}
                self.Melody.Score.convolved_labels[self.Melody.ID][self.ID] = []
                self.Melody.Score.convolved_labels[self.Melody.ID][self.ID].append(Label(parameters,
                     reference = self.processedData.convolved[-1], Phrase=self, 
                     Spectrum=self.processedData.convolved[-1], Melody=self.Melody, 
                     Score=self.Melody.Score))


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
        if not(label == None):
            spectrum = label.Spectrum
            hdr = spectrum.header.copy()
            spectrum.preserve(continuum=False)
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
    def __init__(self, phrases = [], filename=None, name=None, header=None, Score=None):
        self.phrases = phrases
        self.header = header
        self.nPhrases = len(self.phrases)
        self.filename = filename
        self.selectedPhrases = [False for i in range(self.nPhrases)]
        self.muted = True
        self.name = name
        self.Score = Score   # Pointer to Score which contains this melody (optional)
        self.ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        
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

    def record(self, labels=None, outdir='', basename = None):
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
                    filename = outdir+basename+"_T%d_G%.2f_B%.2f_R%d_V%.2f.fits" % (
                        self.Teff, self.logg, self.B, R, vsini)
                print("Recording record \'%s\' to disk" % label.parameters["NAME"])
                label.Phrase.saveConvolved(label=label, filename=filename, header=self.header,
                        wlStart = label.parameters["WLSTART"],
                        wlStop = label.parameters["WLSTOP"])
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
    def __init__(self, phrases = [], filename=None, name=None, header=None, Score=None):
        super(ObservedMelody, self).__init__(phrases=phrases, filename=filename, 
                name=name, header=header, Score=Score)
        if len(self.phrases) == 0:
            if filename != None:
                self.loadMelody()
            else:
                self.phrases = []
                self.nPhrases = 0


    @classmethod
    def fromFile(self, filename=None, name=None, Score=None):
        info = pyfits.info(filename, output='')
        nPhrases = len(info)-1
        header = pyfits.getheader(filename, ext=0)
        phrases = []
        for i in range(nPhrases):
            hdr = pyfits.getheader(filename, ext=i+1)
            phrases.append(ObservedPhrase.fromFile(header=hdr, filename=filename, ext=i+1))

        return self(phrases=phrases, filename=filename, name=name, header=header, Score=Score)
        
    def loadMelody(self):
        info = pyfits.info(self.filename, output='')
        nSpectra = len(info)-1
        self.header = pyfits.getheader(self.filename, ext=0, memmap=False)
        self.phrases = []
        self.contents = self.header.get("SPECTRUM_CONTENTS")
        if self.contents == None:
            self.contents = "OBSERVED"
        parameters = {}
        parameters["NAME"] = "OBSERVED SPECTRUM - %s" % self.name
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
                #label = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % (
                #        self.Teff, self.logg, self.B, hdr.get('VSINI'), hdr.get('RESOLVING_POWER'))
                parameters["NAME"] = self.name
                #parameters["VSINI"] = hdr.get('VSINI')
                #parameters["R"] = hdr.get('RESOLVING_POWER')
                #parameters["WLSTART"] = hdr.get('WLSTART')
                #parameters["WLSTOP"] = hdr.get('WLSTOP')
                parameters["SELECTED"] = False
                self.phrases.append(ObservedPhrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1,
                    parameters=parameters, Melody=self))


        self.nPhrases = len(self.phrases)

    def addPhrase(self, phrase = None):
        parameters = {}
        parameters["SELECTED"] = True   #TODO: What is done with this parameters object?  Nothing!
        self.phrases.append(phrase)
        self.nPhrases = len(self.phrases)

    def loadData(self):
        for phrase in self.phrases:
            phrase.loadData()
            
    def rehearse(self, **kwargs):
        return

    def selectPhrases(self, wlRange=[], exact=False, selectAll=False):
        for phrase in self.phrases:
            if selectAll:
                for observed in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = True
            else:
                for convolved in phrase.convolvedData:
                    convolved.label.parameters["SELECTED"] = phrase.inWlRange(wlStart=wlRange[0],
                                     wlStop=wlRange[1], exact=exact)

    def perform(self):
        spectra = []
        for i in range(self.nPhrases):
            spectra.append(self.phrases[i].listen())
            spectra.append(sp)

        return spectra

    def record(self, labels=[], basename = None):
        """
        Moog960::Melody::record(self, labels=[], basename=None)
        
        record() saves the selected phrases of a melody to disk
        """
        for label in labels:
            print("Recording record \'%s\' to disk" % label.parameters["NAME"])
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
    def __init__(self, phrases = [], filename=None, name=None, header=None, Score=None):
        super(SyntheticMelody, self).__init__(phrases = phrases, filename=filename,
                name=name, header=header, Score=Score)
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
        try:
            parameters["NAME"] = "T = %dK log g = %.1f B = %.2f kG" % (self.Teff, self.logg, self.B)
            parameters["TEFF"] = self.Teff
            parameters["LOGG"] = self.logg
            parameters["BFIELD"] = self.B
        except:
            parameters["NAME"] = "Synthetic Spectrum"
            parameters["TEFF"] = 0
            parameters["LOGG"] = 0
            parameters["BFIELD"] = 0
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
                    name = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI'))
                    parameters["NAME"] = name
                    parameters["VSINI"] = hdr.get('VSINI')
                    parameters["WLSTART"] = hdr.get('WLSTART')
                    parameters["WLSTOP"] = hdr.get('WLSTOP')
                    parameters["SELECTED"] = False
                if self.contents == 'CONVOLVED':
                    name = "T = %dK log g = %.1f B = %.2f kG vsini = %.2f km/s R = %d" % (
                            self.Teff, self.logg, self.B, hdr.get('VSINI'), hdr.get('RESOLVING_POWER'))
                    parameters["NAME"] = name
                    parameters["VSINI"] = hdr.get('VSINI')
                    parameters["R"] = hdr.get('RESOLVING_POWER')
                    parameters["WLSTART"] = hdr.get('WLSTART')
                    parameters["WLSTOP"] = hdr.get('WLSTOP')
                    parameters["SELECTED"] = False
                self.phrases.append(SyntheticPhrase.fromFile(hdr, data=None,
                    filename=self.filename, ext=i+1, diskInt='BEACHBALL',
                    sourceType=self.contents, parameters=parameters, Melody=self))


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

    def rehearse(self, vsini = 0.0, R = 0, observedWl = None):
        convolved = []
        for phrase in self.Score.raw_labels[self.ID].keys():
            self.Score.raw_labels[self.ID][phrase][0].Phrase.rehearse(vsini=vsini, R=R,
                        observedWl=observedWl)

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
            newPhrase = SyntheticPhrase(convolvedData=[mergedSpectra], Melody=self)
            mergedLabel.addReference(Phrase=newPhrase)
            self.Score.appendLabel(mergedLabel)
            self.addPhrases(phrases = [newPhrase])
        
        return mergedLabel
        

class Score( object ):
    """
        This Score object contains many melodies.
    """
    def __init__(self, melodies = [], directory='', observed=None, suffix='raw', files=[]):
        self.ID = ''.join(random.choice(string.ascii_letters) for _ in range(10))
        self.syntheticMelodies = melodies
        self.directory = directory
        self.observed = observed
        self.suffix = suffix
        self.files = files
        self.loadMelodies()
        self.computeGridPoints()
        #self.getMelodyParams(retLabels = False)

    def loadMelodies(self):
        if len(self.files) == 0:
            melodyFiles = glob.glob(self.directory+'*'+self.suffix+'.fits')
        else:
            melodyFiles = self.files
        self.syntheticMelodies = []
        self.raw_labels = {}
        self.interpolated_labels = {} 
        self.integrated_labels = {}
        self.convolved_labels = {}
        self.observed_labels = {}
        for melody in melodyFiles:
            print("%s" % melody)
            self.syntheticMelodies.append(SyntheticMelody(filename=melody, Score=self))

        if not(self.observed==None):
            self.ObservedMelodies = [ObservedMelody(filename=self.observed, 
                Score=self)]
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
                    "NAME", "R"]
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
        for raw in self.getLabels(keySignature='RAW'):
            for key in keys:
                params[key].append(raw.parameters[key])
        for key in params.keys():
            self.rawGridPoints[key] = numpy.unique(params[key])

        # interpolated GridPoints
        params = {}
        self.interpolatedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD"]
        for key in keys:
            params[key] = []
        for interpolated in self.getLabels(keySignature="INTERPOLATED"):
            for key in keys:
                params[key].append(interpolated.parameters[key])
        for key in params.keys():
            self.interpolatedGridPoints[key] = numpy.unique(params[key])

        # integrated GridPoints
        params = {}
        self.integratedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD"]
        for key in keys:
            params[key] = []
        for integrated in self.getLabels(keySignature="INTEGRATED"):
            for key in keys:
                params[key].append(integrated.parameters[key])
        for key in params.keys():
            self.integratedGridPoints[key] = numpy.unique(params[key])

        # convolved GridPoints
        params = {}
        self.convolvedGridPoints = {}
        keys = ["TEFF", "LOGG", "BFIELD", "WLSTART", "WLSTOP", "SELECTED", "VSINI", 
                "NAME", "R"]
        for key in keys:
            params[key] = []
        #for convolved in self.convolved_labels:
        for convolved in self.getLabels(keySignature='CONVOLVED'):
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
        for observed in self.getLabels(keySignature='OBSERVED'):
            for key in keys:#observed.parameters.keys():
                params[key].append(observed.parameters[key])
        for key in params.keys():
            self.observedGridPoints[key] = numpy.unique(params[key])

    def getLabels(self, keySignature='CONVOLVED', selected=False):
        retval = []
        if keySignature == 'RAW':
            for melody in self.raw_labels.keys():
                for phrase in self.raw_labels[melody].keys():
                    for raw in self.raw_labels[melody][phrase]:
                        if selected:
                            if raw.parameters["SELECTED"]:
                                retval.append(raw)
                        else:
                            retval.append(raw)

        if keySignature == 'INTERPOLATED':
            for melody in self.interpolated_labels.keys():
                for phrase in self.interpolated_labels[melody].keys():
                    for interpolated in self.interpolated_labels[melody][phrase]:
                        if selected:
                            if interpolated.parameters["SELECTED"]:
                                retval.append(interpolated)
                        else:
                            retval.append(interpolated)

        if keySignature == 'INTEGRATED':
            for melody in self.integrated_labels.keys():
                for phrase in self.integrated_labels[melody].keys():
                    for integrated in self.integrated_labels[melody][phrase]:
                        if selected:
                            if integrated.parameters["SELECTED"]:
                                retval.append(integrated)
                        else:
                            retval.append(integrated)

        if keySignature == 'CONVOLVED':
            for melody in self.convolved_labels.keys():
                for phrase in self.convolved_labels[melody].keys():
                    for convol in self.convolved_labels[melody][phrase]:
                        if selected:
                            if convol.parameters["SELECTED"]:
                                retval.append(convol)
                        else:
                            retval.append(convol)

        if keySignature == 'OBSERVED':
            for melody in self.observed_labels.keys():
                for phrase in self.observed_labels[melody].keys():
                    for observed in self.observed_labels[melody][phrase]:
                        if selected:
                            if observed.parameters["SELECTED"]:
                                retval.append(observed)
                        else:
                            retval.append(observed)

        return retval

    def appendLabel(self, label=None, keySignature='CONVOLVED'):
        if keySignature == 'RAW':
            self.raw_labels[label.parameters["MELODY"]][label.parameters["PHRASE"]].append(label)
        if keySignature == 'INTEGRATED':
            self.integrated_labels[label.parameters["MELODY"]][label.parameters["PHRASE"]].append(label)
        if keySignature == 'INTERPOLATED':
            self.interpolated_labels[label.parameters["MELODY"]][label.parameters["PHRASE"]].append(label)
        if keySignature == 'CONVOLVED':
            self.convolved_labels[label.parameters["MELODY"]][label.parameters["PHRASE"]].append(label)
        if keySignature == 'OBSERVED':
            self.observed_labels[label.parameters["MELODY"]][label.parameters["PHRASE"]].append(label)

    def selectMelodies(self, wlRange=[], exact=False, keySignature='CONVOLVED'):
        if keySignature=='RAW':
            for melody in self.raw_labels.keys():
                for labels in self.raw_labels[melody].keys():
                    for label in self.raw_labels[melody][labels]:
                        label.parameters["SELECTED"] = label.Phrase.inWlRange(wlStart=wlRange[0], 
                                                       wlStop=wlRange[1], exact=exact)

        if keySignature=='INTERPOLATED':
            for melody in self.interpolated_labels.keys():
                for labels in self.interpolated_labels[melody].keys():
                    for label in self.interpolated_labels[melody][labels]:
                        label.parameters["SELECTED"] = label.Phrase.inWlRange(wlStart=wlRange[0], 
                                                       wlStop=wlRange[1], exact=exact)

        if keySignature=='INTEGRATED':
            for melody in self.integrated_labels.keys():
                for labels in self.integrated_labels[melody].keys():
                    for label in self.integrated_labels[melody][labels]:
                        label.parameters["SELECTED"] = label.Phrase.inWlRange(wlStart=wlRange[0], 
                                                       wlStop=wlRange[1], exact=exact)

        if keySignature=='CONVOLVED':
            for melody in self.convolved_labels.keys():
                for labels in self.convolved_labels[melody].keys():
                    for label in self.convolved_labels[melody][labels]:
                        label.parameters["SELECTED"] = label.Phrase.inWlRange(wlStart=wlRange[0], 
                                                       wlStop=wlRange[1], exact=exact)

        if keySignature=='OBSERVEDED':
            for melody in self.observed_labels.keys():
                for labels in self.observed_labels[melody].keys():
                    for label in self.observed_labels[melody][labels]:
                        label.parameters["SELECTED"] = label.Phrase.inWlRange(wlStart=wlRange[0], 
                                                       wlStop=wlRange[1], exact=exact)

    def selectEnsemble(self, selectedLabels=[], keySignature='CONVOLVED'):
        if keySignature=="CONVOLVED":
            convolved = self.getLabels(keySignature)
            for label in convolved:
                label.parameters["SELECTED"] = False
        for selected in selectedLabels:
            selected.parameters["SELECTED"] = True

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
        parameters = []
        if len(selectedLabels) == 0:
            for melody in self.syntheticMelodies:
                spectrum = []
                params = []
                if keySignature == "CONVOLVED":
                    for convolved in self.getLabels(keySignature='CONVOLVED'):
                        if convolved.parameters["SELECTED"]:
                            spectrum.append(convolved.Spectrum)
                            params.append(convolved)
                if len(spectrum) > 0:
                    spectra.append(spectrum)
                    parameters.append(params)
        else:
            for label in selectedLabels:
                if label.parameters["SELECTED"]:
                     spectra.append(label.Spectrum)
                     parameters.append(label)

        return spectra, parameters


    def listen(self, keySignature='OBSERVED'):   # gets the observed spectrum
        self.compositeObserved = []
        spectra = []
        labels = []
        for melody in self.observed_labels.keys():
            for phrase in self.observed_labels[melody].keys():
                for label in self.observed_labels[melody][phrase]:
                    spectra.append(label.Phrase.listen())
                    labels.append(label)

        spectra = numpy.array(spectra)
        labels = numpy.array(labels)

        order = numpy.argsort(labels)
        spectra = spectra[order]
        labels = labels[order]

        mergedSpectra = spectra[0][0]
        mergedLabel = labels[0]
        for sp, l in zip(spectra[1:], labels[1:]):
            mergedSpectra = mergedSpectra.mergeSpectra(sp[0])
            mergedLabel = mergedLabel.merge(l)

        mergedLabel.parameters["SELECTED"] = True
        mergedLabel.addReference(Spectrum = mergedSpectra)
        if keySignature=='OBSERVED':
            newPhrase = ObservedPhrase(observedData=[mergedSpectra], Melody=self)
            mergedLabel.addReference(Phrase=newPhrase)
            mergedLabel.addReference(Score=self)
            self.appendLabel(mergedLabel, keySignature=keySignature)

        self.compositeObservedLabel = mergedLabel

        return mergedSpectra, mergedLabel

    def record(self, selected=[], filename='', keySignature='CONVOLVED', dI = False, 
              primaryHeaderKWs=None):
        if len(selected) > 0:
            for label in selected:
               label.Phrase.record(filename=filename)
        else:
            if keySignature == 'OBSERVED':
                for label in self.getLabels(keySignature=keySignature):
                    label.Phrase.record(filename=filename, dI=dI, primaryHeaderKWs=primaryHeaderKWs)
            if keySignature == 'CONVOLVED':
                for label in self.getLabels(keySignature=keySignature):
                    label.Phrase.record(filename=filename)
            

    def master(self, selectedLabels=[], keySignature='CONVOLVED'):
        """
           Moog960::SyntheticMelody::master()
           This routine merges all phrases into a single master phrase
        """
        mastered = []
        for synthetic in self.syntheticMelodies:
            mastered.append(synthetic.master(keySignature=keySignature))
            mastered[-1].parameters["SELECTED"] = True

        return mastered


    def blend(self, desiredParameters={}, appendTheBlend=True):
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
        for label in self.getLabels(keySignature='CONVOLVED'):
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
            p.addReference(Spectrum=sp, Score=self)
            phrases.append(SyntheticPhrase(convolvedData=[sp], diskInt='BEACHBALL'))
            #self.appendLabel(label=p, keySignature='CONVOLVED')
        header=pyfits.Header()
        for key in desiredParameters.keys():
            header.set(key, desiredParameters[key])
        header.set("SPECTRUM_CONTENTS", "CONVOLVED")
        blended = SyntheticMelody(phrases=phrases, header=header, Score=self)
        for p in interpolatedParams:
            p.addReference(Melody=blended)
            self.appendLabel(label=p, keySignature='CONVOLVED')

        if appendTheBlend:
            self.syntheticMelodies.append(blended)
        return blended, interpolatedParams
        
    def calc_lnlike(self, Teff=0.0, logg=0.0, Bfield=0.0, rv=0.0, ax=None):
        desiredParams = {}
        desiredParams["TEFF"] = Teff
        desiredParams["LOGG"] = logg
        desiredParams["BFIELD"] = Bfield
        print rv
        blended, blendedLabels = self.blend(desiredParameters=desiredParams,
                  appendTheBlend=False)
        print "Blend Finished!"
        blendedLabels[0].Spectrum.rv(rv)
        print "RV Finished!"
        blendedLabels[0].Spectrum.bin(self.compositeObservedLabel.Spectrum.wl, pad=0.0)
        print "Binning Finished!"
        difference = blendedLabels[0].Spectrum - self.compositeObservedLabel.Spectrum
        print "Difference Finished!"
        lnlike = -0.5*numpy.sum( 
                  (difference/self.compositeObservedLabel.Spectrum).flux_I**2.0)
        print "ln_likelihood Finished!"
        print lnlike
        if ax != None:
            ax.clear()
            blendedLabels[0].Spectrum.plot(ax=ax)
            self.compositeObservedLabel.Spectrum.plot(ax=ax)
            ax.figure.show()
            raw_input()
        return lnlike

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
