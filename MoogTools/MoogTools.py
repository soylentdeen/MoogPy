import scipy
import scipy.interpolate
import numpy
import os
import sys
import SpectralTools
import AstroUtils
import astropy.io.fits as pyfits
import time
import Moog960

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
            self.rho[i] = self.pgas[i]*self.molweight[i]*1.6606e-24/(1.38054e-16*self.  T[i])
            self.kapref[i] = self.kaprefmass[i] * self.rho[i]
        self.mt = float(data.readline().split()[0])
        data.close()

class progressBar( object ):
    def __init__(self, start, stop):
        self.start_value = start
        self.stop_value = stop

    def start(self):
        self.currentValue = self.start_value
        self.numBlocks = 0
        print("Starting Synthesis!")

    def update(self, value):
        self.currentValue = value
        self.percentComplete = (self.currentValue-self.start_value)/(self.stop_value
                -self.start_value)*100
        if int(self.percentComplete/5) > self.numBlocks:
            sys.stdout.write("\r%.2f [%s%s] %.2f - %.2f" % (self.start_value,
                '#'*int(self.percentComplete/5),' '*(20-int(self.percentComplete/5)),
                self.stop_value, self.currentValue))
            sys.stdout.flush()
            self.numBlocks = int(self.percentComplete/5)


class periodicTable( object ):
    def __init__(self):
        self.Zsymbol_table = {}
        self.DissE_table = {}
        df = open(os.path.dirname(os.path.realpath(__file__))+'/MOOGConstants.dat', 'r')
        for line in df.readlines():
            l = line.split('-')
            self.Zsymbol_table[int(l[0])] = l[1].strip()
            self.Zsymbol_table[l[1].strip()] = int(l[0])
            if len(l) > 3:
                self.DissE_table[int(l[0])] = float(l[3])
                self.DissE_table[l[1].strip()] = float(l[3])
        df.close()

    def translate(self, ID):
        retval = self.Zsymbol_table[ID]
        return retval

    def DissE(self, ID):
        retval = self.DissE_table[ID]
        return retval

class HITRAN_Dictionary( object ):
    def __init__(self):

        # HITRAN codes:
        #  5 : CO
        # 13 : OH
        self.isotopes = {5:{1:608.01216,2:608.01316,3:608.01218,4:608.01217,5:608.01318,6:608.01317},
                13:{1:108.01160,2:108.01180,3:108.02160}}
        self.DissE = {5:11.10, 13:4.412}

class MoogStokes( object ):
    def __init__(self, configurationFile, fileBase="moogstokes", **kwargs):
        """
        MoogStokes:  A python wrapper for MoogStokes

        Input:
            configurationFile : The name of a file containing configuration parameters
            fileBase : the base name for the MoogStokes parameter file/linelists

        Member variables
            config - a dictionary created from lines contained in the configuration file
            lineList - a LineList object containing lines within the defined wavelength range
            parameterFile - a ParameterFile object useful for writing Moog-Readable parameter files
            T - Temperature in Kelvin
            logg - Surface gravity
            B - Magnetic field strength in kG

            fileName - name of the MoogStokes parameter file
            fileBase - base of the filename
            wave = array containing wavelengths of the spectrum
            flux = array containing the flux of the composite emergent spectrum
        """
        self.config = AstroUtils.parse_config(configurationFile)
        if "MODELFILE" in kwargs.keys():
            self.modelFile = kwargs["MODELFILE"]
            self.T = 0.0
            self.logg = 0.0
            self.B = 0.0
        elif "modelFile" in self.config.keys():
            self.modelFile = self.config["modelFile"]
            self.T = 0.0
            self.logg = 0.0
            self.B = 0.0
        else:
            self.T = self.config["Teff"]
            self.logg = self.config["logg"]
            self.B = self.config["Bfield"]
            self.modelFile = None
        
        if "diskInt" in kwargs.keys():
            self.diskInt = kwargs["diskInt"]
        elif "diskInt" in self.config.keys():
            self.diskInt = self.config["diskInt"]
        else:
            self.diskInt = "BEACHBALL"

        if self.diskInt == "TENNISBALL":
            self.LineListFormat="MOOGSCALAR"
        else:
            self.LineListFormat = "MOOGSTOKES"

        if "MoogSandbox" in self.config.keys():
            self.MoogSandbox = self.config["MoogSandbox"]
        else:
            self.MoogSandbox = ''

        if "wlStart" in kwargs.keys():
            self.config["wlStart"] = kwargs["wlStart"]
        if "wlStop" in kwargs.keys():
            self.config["wlStop"] = kwargs["wlStop"]
        self.fileName = self.MoogSandbox+fileBase+'.par'
        self.fileBase = fileBase
        self.lineList = LineList(self, self.config)
        self.parameterFile = ParameterFile(self, self.config)
        self.Spectra = []
        self.logtau = []

        if "moogInstance" in kwargs.keys():
            MoogInstance = kwargs["moogInstance"].upper()
            if MoogInstance == "ALPHA":
                import MoogStokesPy_Alpha
                self.MoogPy = MoogStokesPy_Alpha
            elif MoogInstance == "BRAVO":
                import MoogStokesPy_Bravo
                self.MoogPy = MoogStokesPy_Bravo
            elif MoogInstance == "CHARLIE":
                import MoogStokesPy_Charlie
                self.MoogPy = MoogStokesPy_Charlie
            elif MoogInstance == "DELTA":
                import MoogStokesPy_Delta
                self.MoogPy = MoogStokesPy_Delta
        else:
            import MoogStokesPy
            self.MoogPy = MoogStokesPy
        self.MoogPy.charstuff.moogpath = '%-60s'%os.environ.get('MOOGSTOKESSOURCE')
        self.MoogPy.recorder = self.recorder
        self.MoogPy.stokesrecorder = self.stokesrecorder
        self.MoogPy.beachball = self.beachball
        self.MoogPy.diskoball = self.diskoball
        self.MoogPy.fluxtracer = self.fluxtracer
        self.MoogPy.tennisball = self.tennisball

        if "progressBar" in kwargs.keys():
            if kwargs["progressBar"] == True:
                self.progressBar = progressBar(self.config["wlStart"], self.config["wlStop"])
            else: self.progressBar = None
        else:
            self.progressBar = None

    def fluxtracer(self, logtau, dtau, Stokes, continuum):
        self.logtau.append(logtau)
        self.flux_I.append(Stokes[0])
        self.flux_Q.append(Stokes[1])
        self.flux_U.append(Stokes[2])
        self.flux_V.append(Stokes[3])
        self.continuum.append(continuum)

    def recorder(self, x, y):
        self.Spectra[-1].wl.append(x)
        self.Spectra[-1].flux_I.append(1.0-y)
    
    def stokesrecorder(self, i, wave, Stokes, continuum):
        index = i-1
        self.Spectra[index].wl.append(wave)
        self.Spectra[index].flux_I.append(Stokes[0])
        self.Spectra[index].flux_Q.append(Stokes[1])
        self.Spectra[index].flux_U.append(Stokes[2])
        self.Spectra[index].flux_V.append(Stokes[3])
        self.Spectra[index].continuum.append(continuum)
        if i == 1:
            if self.progressBar != None:
                self.progressBar.update(wave)

    def prepareFluxes(self):
        for i in range(self.ncells):
            header = pyfits.Header()
            header.set('CREATION_TIME', time.ctime())
            try:
               header.set('CREATION_USER', os.getlogin())
            except:
               header.set('CREATION_USER', os.getenv('USER'))
            header.set('CREATION_MACHINE', os.uname()[1])
            header.set('MOOGVERSION', self.MoogPy.charstuff.moogversion.tostring())
            header.set('WLSTART', self.config["wlStart"])
            header.set('WLSTOP', self.config["wlStop"])
            header.set('CELL', self.cells[i])
            if self.diskInt == 'TENNISBALL':
                header.set('SPECTRUM_TYPE', 'MOOG disk-integrated Spectrum')
                self.Spectra.append(SpectralTools.Spectrum(wl = [],
                    I = [], header=header,
                    spectrum_type='MOOG DISK INTEGRATED'))
            elif self.diskInt == 'BEACHBALL':
                header.set('PHI_ANGLE', float('%.4f'% self.phi_angle[i]))
                header.set('MU', float('%.4f'%self.mus[i]))
                header.set('DISKFLAG', self.diskflag)
                header.set('SPECTRUM_TYPE', 'BeachBall Emergent Spectrum')
                self.Spectra.append(SpectralTools.Spectrum(wl = [],
                    I = [], Q = [], U = [], V = [], continuum = [], header=header, 
                    spectrum_type='MOOG EMERGENT'))
            elif self.diskInt == 'DISKOBALL':
                header.set('PHI_ANGLE', float('%.4f'% self.phi_angle[i]))
                header.set('CHI_ANGLE', float('%.4f'% self.chi_angle[i]))
                header.set('AZIMUTH', float('%.4f' % self.azimuth[i]))
                header.set('LONGITUDE', float('%.4f' % self.longitude[i]))
                header.set('NRINGS', float('%d' % self.nrings))
                header.set('INCLINATION', float('%.4f' % self.inclination))
                header.set('POSANGLE', float('%.4f' % self.position_angle))
                header.set('MU', float('%.4f'%self.mus[i]))
                header.set('DISKFLAG', self.diskflag)
                header.set('SPECTRUM_TYPE', 'DiskoBall Emergent Spectrum')
                self.Spectra.append(SpectralTools.Spectrum(wl = [],
                    I = [], Q = [], U = [], V = [], continuum = [], header=header, 
                    spectrum_type='MOOG EMERGENT'))

    def tennisball(self):
        self.ncells = 1
        self.cells = numpy.arange(1)
        self.LineListFormat = "MOOGSCALAR"
        self.prepareFluxes()

    def beachball(self):
        self.diskflag = 1
        self.ncells = 7
        try:
            self.deltav = self.config["deltav"]
        except:
            self.deltav = 0.1
        self.cells = numpy.arange(7)
        self.phi_angle = self.MoogPy.angles.phi_angle.copy()
        self.mus = self.MoogPy.angles.mus.copy()
        self.prepareFluxes()

    def diskoball(self):
        self.diskflag = 0
        self.ncells = self.MoogPy.angles.ncells
        self.nrings = self.MoogPy.angles.nrings
        self.inclination = self.MoogPy.angles.inclination
        self.position_angle = self.MoogPy.angles.position_angle
        self.phi_angle = self.MoogPy.angles.phi_angle.copy()
        self.chi_angle = self.MoogPy.angles.chi_angle.copy()
        self.azimuth = self.MoogPy.angles.azimuth.copy()
        self.longitude = self.MoogPy.angles.longitude.copy()
        self.mus = self.MoogPy.angles.mus.copy()
        self.prepareFluxes()

    def finishSpectra(self):
        for spectrum in self.Spectra:
            spectrum.preserve(I=spectrum.flux_I!=None, Q=spectrum.flux_Q!=None,
                    U=spectrum.flux_U!=None, V=spectrum.flux_V!=None, 
                    continuum=spectrum.continuum!=None)

    def run(self, saveRaw=False, **kwargs):
        self.lineList.setBfield(self.B)
        self.lineList.writeLineLists(parent=self, mode=self.LineListFormat, **kwargs)
        self.parameterFile.setName(self.fileBase)
        self.parameterFile.setModel(teff=self.T, logg=self.logg, 
                modelFile=self.modelFile)
        self.parameterFile.writeParFile()
        self.MoogPy.charstuff.fparam = self.fileName.ljust(80)
        self.MoogPy.atmos.linecount = 0
        if self.progressBar != None:
            self.progressBar.start()
        self.MoogPy.moogstokessilent()
        self.finishSpectra()
        self.Phrase = Moog960.SyntheticPhrase(rawData=self.Spectra,
                diskInt=self.diskInt)
        if saveRaw:
            filename = self.config["outdir"]+self.config["outbase"]+'_T%d_G%.2f_B%.2f_raw.fits'%(self.config["Teff"],
                    self.config["logg"], self.config["Bfield"])
            PHKWs = {"BFIELD":self.config["Bfield"], "TEFF":self.config["Teff"], "LOGG":self.config["logg"]}
            self.Phrase.saveRaw(filename=filename, primaryHeaderKWs=PHKWs)

    def trace(self, save=False):
        self.lineList.setBfield(self.B)
        self.lineList.writeLineLists(parent=self, mode="MOOGSTOKES")
        self.MoogPy.charstuff.fparam = self.fileName.ljust(80)
        self.parameterFile.setName(self.fileBase)
        self.parameterFile.setModel(teff = self.T, logg = self.logg)
        self.parameterFile.writeParFile()
        self.MoogPy.atmos.linecount = 0
        self.MoogPy.moogstokessilent()
        if save:
            out = pyfits.PrimaryHDU([self.logtau, self.flux_I, self.flux_Q, self.flux_U, self.flux_V, self.continuum])
            out.header.set('WAVE', self.config["wlProbe"])
            out.header.set('BFIELD', self.config["Bfield"])
            out.header.set('TEFF', self.config["Teff"])
            out.header.set('LOGG', self.config["logg"])
            out.writeto(self.config["outdir"]+self.config["outbase"]+'_'+str(self.config["wlProbe"])+'.fits', clobber=True)
        return self.logtau, self.flux_I, self.flux_Q, self.flux_U, self.flux_V, self.continuum

class ParameterFile( object ):
    def __init__(self, parent, config, **kwargs):
        self.parent = parent
        self.config = config
        self.moogParCfgFile = config['moog_Parameters']
        self.moogPars = AstroUtils.parse_config(self.moogParCfgFile)
        self.synlimits = numpy.array([config['wlStart'], config['wlStop'], 
                self.moogPars['synlimits_c'], self.moogPars['synlimits_d']])
        if "wlProbe" in self.config.keys():
            self.wlProbe = self.config["wlProbe"]
        else:
            self.wlProbe = None
        if "PARFILENAME" in kwargs:
            self.parFileName = kwargs["PARFILENAME"]
        else:
            self.parFileName = self.moogPars['parFileName']
        if "atmos_dir" in self.moogPars.keys():
            atmos_dir = self.moogPars["atmos_dir"]
            self.moogPars["atmos_dir"] = os.environ.get('MOOGPYDATAPATH')+atmos_dir
        self.mode = self.moogPars['mode']
        self.labels = {'terminal':'x11',
                      'strong':1, 
                      'atmosphere':1, 
                      'molecules':2,
                      'lines':0,
                      'damping':1,
                      'freeform':2,
                      'flux/int':0,
                      'diskflag':1,
                      'testflag':0}
        self.file_labels = {'summary_out':'./Output/summary.out',
                            'standard_out':'./Output/out1',
                            'smoothed_out':'./Output/smoothed.out',
                            'model_in':'',
                            'lines_in':config['Weak_FileName'],
                            'stronglines_in':config['Strong_FileName']}
                            #'out_dir':'',
        for l in self.labels:
            if l in self.moogPars:
                self.labels[l] = self.moogPars[l]
            
        for fl in self.file_labels:
            if fl in self.moogPars:
                self.file_labels[fl] = self.moogPars[fl]

    def setName(self, name):
        self.file_labels['lines_in'] = self.parent.MoogSandbox + self.config['Weak_FileName']+'_'+name
        self.file_labels['stronglines_in'] = self.parent.MoogSandbox + self.config['Strong_FileName']+'_'+name
        self.parFileName = self.parent.MoogSandbox + name+'.par'
        
    def setModel(self, teff=0.0, logg=0.0, modelFile=None):
        if modelFile==None:
            self.file_labels["model_in"] = os.environ.get('MOOGPYDATAPATH')+ \
                    'Atmospheres/MARCS/MARCS_T'+ str(int(teff))+'_G'+ \
                    str(logg)+'_M0.0_t1.0.md'
        else:
            self.file_labels["model_in"] = os.environ.get('MOOGPYDATAPATH') + \
                    'Atmospheres/' + modelFile

    def writeParFile(self):
        pf = open(self.parFileName, 'w')

        pf.write(self.mode+'\n')
        for fl in self.file_labels:
            pf.write(fl+'   \''+self.file_labels[fl]+'\'\n')
        for l in self.labels:
            pf.write(l+'    '+str(self.labels[l])+'\n')

        if self.wlProbe:
            pf.write('dipstick  %.3f\n' % self.wlProbe)

        pf.write('synlimits\n')
        pf.write('               %.2f %.2f %.3f %.2f\n' % 
                (self.synlimits[0], self.synlimits[1], 
                self.synlimits[2], self.synlimits[3]))
        pf.close()
        


class LineList( object ):
    def __init__(self, parent, config):
        # Load in configuration file
        self.parent = parent
        self.MoogPyDataPath = os.environ.get('MOOGPYDATAPATH')
        self.strong_file = self.MoogPyDataPath+config['strong_file']
        self.VALD_list = self.MoogPyDataPath+config['VALD_file']
        self.gf_corrections = self.MoogPyDataPath+config['gf_file']
        self.wlStart = config['wlStart']
        self.wlStop = config['wlStop']
        self.Bfield = config['Bfield']/10.0
        self.sfn = parent.MoogSandbox + config['Strong_FileName']
        self.wfn = parent.MoogSandbox + config['Weak_FileName']
        self.applyCorrections = config['applyCorrections']

        self.readInLineLists()
        self.nStrong = len(self.strongLines)
        self.numLines = self.nStrong+len(self.weakLines)
        self.dummyLine = dummy_Line()

    def setBfield(self, B):
        self.Bfield = B/10.0
        for i in range(self.nStrong):
            self.strongLines[i].zeeman_splitting(B=self.Bfield)
        for i in range(self.numLines - self.nStrong):
            if not(self.weakLines[i].DissE):
                self.weakLines[i].zeeman_splitting(B=self.Bfield)

    def readInLineLists(self):
        self.parse_new_VALD()

    def parse_new_VALD(self):
        pt = periodicTable()
    
        if self.applyCorrections:
            self.corrected = []
            for line in open(self.gf_corrections, 'r'):
                self.corrected.append(MOOG_Line(line))
    
        strong = []
        for line in open(self.strong_file, 'r'):
            l = line.split()
            strong.append([float(l[0]), float(l[1])])
    
        vald_in = open(self.VALD_list, 'r')
        l1 = ''
        l2 = ''
        l3 = ''
        l4 = ''
        self.strongLines = []
        self.weakLines = []
        junk = vald_in.readline()
        junk = vald_in.readline()
        linecounter = 0
        lines = [l1, l2, l3, l4]
        for line in vald_in:
            lines[linecounter] = line
            linecounter += 1
            if linecounter == 4:
                linecounter = 0
                if (lines[0][0] == '\''):
                    wl = float(lines[0].split(',')[1])
                    if ( (wl > self.wlStart) & (wl < self.wlStop) ):
                        current_line = New_VALD_Line(lines, pt)
                        if ( (current_line.expot_lo < 20.0) & 
                                (current_line.species % 1 <= 0.2) ):
                            if self.applyCorrections:
                                for cl in self.corrected:
                                    if ((cl.wl == current_line.wl) & 
                                            (cl.expot_lo == current_line.expot_lo) &
                                            (cl.loggf != -6.0) ):
                                        current_line.loggf = cl.loggf
                                        current_line.zeeman["NOFIELD"][1] = cl.loggf
                                        if ((cl.VdW != 99.0) & (cl.species < 100.0)):
                                            current_line.VdW = cl.VdW
                                            current_line.stark = cl.stark
                                            current_line.radiative = cl.radiative
                            current_line.zeeman_splitting(self.Bfield)
                            species = current_line.species
                            if ( [wl, species] in strong):
                                self.strongLines.append(current_line)
                            else:
                                self.weakLines.append(current_line)
    
    def getWl(self, index):
        if index < self.nStrong:
            return self.strongLines[index].wl
        else:
            return self.weakLines[index-self.nStrong].wl

    def getGf(self, index, log=False):
        if index < self.nStrong:
            if log:
                return self.strongLines[index].loggf
            else:
                return 10.0**(self.strongLines[index].loggf)
        else:
            if log:
                return self.weakLines[index-self.nStrong].loggf
            else:
                return 10.0**(self.weakLines[index-self.nStrong].loggf)
        
    def getVdW(self, index, log=False):
        if index < self.nStrong:
            if log:
                return self.strongLines[index].VdW
            else:
                return 10.0**(self.strongLines[index].VdW)
        else:
            if log:
                return self.weakLines[index-self.nStrong].VdW
            else:
                return 10.0**(self.weakLines[index-self.nStrong].VdW)

    def getWl(self, index):
        if index < self.nStrong:
            return self.strongLines[index].wl
        else:
            return self.weakLines[index-self.nStrong].wl

    def perturbGf(self, index, delta, push=False):
        if index < self.nStrong:
            self.strongLines[index].modifyGf(delta, push=push)
        else:
            self.weakLines[index-self.nStrong].modifyGf(delta, push=push)

    def perturbVdW(self, index, delta, push=False):
        if index < self.nStrong:
            self.strongLines[index].modifyVdW(delta, push=push)
        else:
            self.weakLines[index-self.nStrong].modifyVdW(delta, push=push)

    def saveLineList(self, mode="MOOGSCALAR", filename='', changed=False):
        outfile = open(filename, 'w')
        for strongLine in self.strongLines:
            strongLine.dump(out=outfile, mode=mode, changed=changed)
        for weakLine in self.weakLines:
            weakLine.dump(out=outfile, mode=mode, changed=changed)
        outfile.close()
        self.sort_file(filename)

    def writeLineLists(self, lineIndex=-1, partial=False, mode="MOOGSCALAR", parent=None):
        if parent:
            strongFile = self.sfn + '_' + parent.fileBase
            weakFile = self.wfn + '_' + parent.fileBase
            moogPointer = parent.MoogPy
        else:
            strongFile = self.sfn
            weakFile = self.wfn
            moogPointer = self.parent.MoogPy
        self.parent.MoogPy.linex.start = self.wlStart
        self.parent.MoogPy.linex.sstop = self.wlStop
        self.parent.parameterFile.synlimits[0] = self.wlStart
        self.parent.parameterFile.synlimits[1] = self.wlStop
        if lineIndex < 0:     #Normal case, write ALL the lines
            outfile = open(strongFile, 'w')
            for strongLine in self.strongLines:
                strongLine.dump(out=outfile, mode=mode)
            outfile.close()
            self.sort_file(strongFile)

            outfile = open(weakFile, 'w')
            if len(self.weakLines) == 0:
                self.dummyLine.create(self.strongLines[lineIndex].wl, outfile)
            else:
                for weakLine in self.weakLines:
                    weakLine.dump(out=outfile, mode=mode)
            outfile.close()
            self.sort_file(weakFile, weak=True)
            """    This worked previously
            self.parent.parameterFile.writeParFile()
            self.parent.MoogPy.linex.start = self.wlStart
            self.parent.MoogPy.linex.sstop = self.wlStop
            self.parent.MoogPy.linex.nlines = len(self.strongLines)
            self.parent.MoogPy.linex.nstrong = len(self.weakLines)
            """
            moogPointer.linex.start = self.wlStart
            moogPointer.linex.sstop = self.wlStop
            moogPointer.linex.nlines = len(self.strongLines)  #Do I need to count lines or components?
            moogPointer.linex.nstrong = len(self.weakLines)
        elif lineIndex < self.nStrong:   # We want to only print out one line
            #  STRONG LINE
            outfile = open(strongFile, 'w')
            self.strongLines[lineIndex].dump(out=outfile, mode=mode)
            outfile.close()
            self.sort_file(strongFile)

            outfile = open(weakFile, 'w')
            self.dummyLine.create(self.strongLines[lineIndex].wl, outfile)
            outfile.close()
            self.sort_file(weakFile, weak=True)

            # Set Moog varibles
            self.parent.MoogPy.linex.nlines = 1
            self.parent.MoogPy.linex.nstrong = 1
        else:
            # WEAK LINE
            index = lineIndex-self.nStrong
            # Write empty strong line file
            outfile = open(strongFile, 'w')
            outfile.close()

            weakLine = self.weakLines[index]
            wlStart = weakLine.wl - 3.0
            wlStop = weakLine.wl + 3.0
            if partial:
                self.parent.parameterFile.synlimits[0] = wlStart
                self.parent.parameterFile.synlimits[1] = wlStop
                self.parent.MoogPy.linex.start = wlStart
                self.parent.MoogPy.linex.sstop = wlStop
                #self.parent.parameterFile.writeParFile()
                #print self.parent.MoogPy.linex.start
                #print self.parent.parameterFile.synlimits
                #raw_input()
                #self.parent.parameterFile.synlimits[0] = self.wlStart
                #self.parent.parameterFile.synlimits[1] = self.wlStop
            self.parent.MoogPy.linex.nlines = 1
            self.parent.MoogPy.linex.nstrong = 0
            outfile = open(weakFile, 'w')
            if weakLine.loggf > 0:
                self.dummyLine.create(weakLine.wl+0.01, outfile)
            weakLine.dump(out=outfile, mode=mode)
            outfile.close()
            self.sort_file(weakFile, weak=True)

    def sort_file(self, name, weak=False):
        data = open(name, 'r').readlines()
        wl = []
        for line in data:
            wl.append(float(line[0:10]))

        order = numpy.argsort(wl)
        out = open(name, 'w')
        if weak:
            out.write("Weak Lines\n")
        for i in order:
            out.write(data[i])
        out.close()

    def applyCorrection(self, corrections):
        #"""
        for i in range(self.numLines):
            if i < self.nStrong:
                self.strongLines[i].modifyGf(corrections[i], push=True)
                self.strongLines[i].modifyVdW(corrections[i+self.numLines], push=True)
            else:
                self.weakLines[i-self.nStrong].modifyGf(corrections[i], push=True)
                self.weakLines[i-self.nStrong].modifyVdW(corrections[i+self.numLines], push=True)

        self.writeLineLists(parent=self.parent)

    def setLogGfs(self, loggfs):
        for i in range(self.numLines):
            if i < self.nStrong:
                self.strongLines[i].setLogGf(loggfs[i])
            else:
                self.weakLines[i-self.nStrong].setLogGf(loggfs[i])
    
    def setVdWs(self, VdWs):
        for i in range(self.numLines):
            if i < self.nStrong:
                self.strongLines[i].setVdW(VdWs[i])
            else:
                self.weakLines[i-self.nStrong].setVdW(VdWs[i])

    def tossLines(self, indices):
        nTossed = 0
        for index in indices:
            index -= nTossed
            if index < self.nStrong:
                junk = self.strongLines.pop(index)
                self.nStrong -= 1
            else:
                junk = self.weakLines.pop(index-self.nStrong)
            nTossed+=1
            self.numLines -= 1

class Spectral_Line( object ):
    def __init__(self):
        self.wl = None
        self.species = None
        self.expot_lo = None
        self.loggf = None
        self.DissE = None
        self.VdW = None
        self.loggfHistory = []
        self.VdWHistory = []
        self.radiative = None
        self.stark = None
        self.zeeman = {}
        self.transition = None
        self.J_lo = None
        self.J_hi = None
        self.g_lo = None
        self.g_hi = None
        self.g_eff = None
        self.verbose = False
        self.Bfield = 0.0
    
    def __str__(self):
        return "%.1f line at %.4f: log gf = %.3f, VdW = %.3f" % (self.species,
                self.wl, self.loggf, self.VdW)

    def __eq__(self, other):
        return (self.wl == other.wl) & (self.species == 
                other.species) & (self.expot_lo == other.expot_lo)

    def modifyGf(self, delta_loggf, push=False):
        if push:
            self.loggfHistory.append(self.loggf)
        self.loggf = numpy.log10(10.0**self.loggf + delta_loggf)
        if numpy.isnan(self.loggf):
            self.loggf = -6.0
        self.zeeman_splitting()

    def modifyVdW(self, delta_VdW, push=False):
        if not(self.VdW):
            self.VdW = -7.5  # No damping value currently, try a guess
        if push:
            self.VdWHistory.append(self.VdW)
        self.VdW += delta_VdW
        if self.VdW >= -5.0:
            self.VdW = -4.999
        if self.VdW < -9.5:
            self.VdW = -9.5

    def setLogGf(self, loggf):
        self.loggf = loggf
        self.zeeman_splitting()

    def setVdW(self, VdW):
        self.VdW = VdW

    def dump(self, **kwargs):
        if "changed" in kwargs:
            if kwargs["changed"] == True:
                if ((len(self.loggfHistory) == 0) & (len(self.VdWHistory) == 0)):
                    return
                if (self.VdW == -9.5):
                    return
                if (self.loggf == -6.0):
                    return
        if "out" in kwargs:
            out = kwargs["out"]
            if kwargs["mode"].upper() == 'MOOGSTOKES':
                if( (self.expot_lo < 20.0) & (self.species % 1 <= 0.2)):
                    if not(self.DissE):
                        for i in range(len(self.zeeman["PI"][0])):
                            out.write('%10.3f%10s%10.3f%10.5f' %
                               (self.zeeman["PI"][0][i],
                               self.species,self.expot_lo,self.zeeman["PI"][1][i]))
                            if not(self.VdW):
                                out.write('%20s%20.3f'% (' ',0.0))
                            else:
                                out.write('%10.3f%20s%10.3f' %
                                        (self.VdW, ' ', 0.0))
                            if not(self.radiative):
                                out.write('%10.3s'% (' '))
                            else:
                                out.write('%10.3f' %
                                        (self.radiative))
                            if not(self.stark):
                                out.write('%10s\n'% (' '))
                            else:
                                out.write('%10.3f\n' %
                                        (self.stark))
                        for i in range(len(self.zeeman["LCP"][0])):
                            out.write('%10.3f%10s%10.3f%10.5f' %
                               (self.zeeman["LCP"][0][i],
                               self.species,self.expot_lo,self.zeeman["LCP"][1][i]))
                            if not(self.VdW):
                                out.write('%20s%20.3f'% (' ',-1.0))
                            else:
                                out.write('%10.3f%20s%10.3f' %
                                        (self.VdW, ' ', -1.0))
                            if not(self.radiative):
                                out.write('%10.3s'% (' '))
                            else:
                                out.write('%10.3f' %
                                        (self.radiative))
                            if not(self.stark):
                                out.write('%10s\n'% (' '))
                            else:
                                out.write('%10.3f\n' %
                                        (self.stark))
                        for i in range(len(self.zeeman["RCP"][0])):
                            out.write('%10.3f%10s%10.3f%10.5f' %
                               (self.zeeman["RCP"][0][i],
                               self.species,self.expot_lo,self.zeeman["RCP"][1][i]))
                            if not(self.VdW):
                                out.write('%20s%20.3f'% (' ',1.0))
                            else:
                                out.write('%10.3f%20s%10.3f' %
                                        (self.VdW, ' ', 1.0))
                            if not(self.radiative):
                                out.write('%10.3s'% (' '))
                            else:
                                out.write('%10.3f' %
                                        (self.radiative))
                            if not(self.stark):
                                out.write('%10s\n'% (' '))
                            else:
                                out.write('%10.3f\n' %
                                        (self.stark))
                    else:
                        #RCP
                        out.write('%10.3f%10.5f%10.3f%10.3f' %
                                (self.wl, self.species, self.expot_lo,self.loggf))
                        if not(self.VdW):
                            out.write('%10s%10.3f%20.3f' %
                                    (' ',self.DissE, 1.0))
                        else:
                            out.write('%10.3f%10.3f%20.3f' %
                                    (self.VdW, self.DissE, 1.0))
                        if not(self.radiative):
                            out.write('%10.3s'% (' '))
                        else:
                            out.write('%10.3f' %
                                    (self.radiative))
                        if not(self.stark):
                            out.write('%10s\n'% (' '))
                        else:
                            out.write('%10.3f\n' %
                                    (self.stark))
                        #PI
                        out.write('%10.3f%10.5f%10.3f%10.3f' %
                                (self.wl, self.species, self.expot_lo,self.loggf))
                        if not(self.VdW):
                            out.write('%10s%10.3f%20.3f' %
                                    (' ',self.DissE, 0.0))
                        else:
                            out.write('%10.3f%10.3f%20.3f' %
                                    (self.VdW, self.DissE, 0.0))
                        if not(self.radiative):
                            out.write('%10.3s'% (' '))
                        else:
                            out.write('%10.3f' %
                                    (self.radiative))
                        if not(self.stark):
                            out.write('%10s\n'% (' '))
                        else:
                            out.write('%10.3f\n' %
                                    (self.stark))
                        #LCP
                        out.write('%10.3f%10.5f%10.3f%10.3f' %
                                (self.wl, self.species, self.expot_lo,self.loggf))
                        if not(self.VdW):
                            out.write('%10s%10.3f%20.3f' %
                                    (' ',self.DissE, -1.0))
                        else:
                            out.write('%10.3f%10.3f%20.3f' %
                                    (self.VdW, self.DissE, -1.0))
                        if not(self.radiative):
                            out.write('%10.3s'% (' '))
                        else:
                            out.write('%10.3f' %
                                    (self.radiative))
                        if not(self.stark):
                            out.write('%10s\n'% (' '))
                        else:
                            out.write('%10.3f\n' %
                                    (self.stark))
            elif kwargs["mode"].upper() == "MOOGSCALAR":
                if( (self.expot_lo < 20.0) & (self.species % 1 <= 0.2)):
                    if not(self.DissE):
                        out.write('%10.3f%10s%10.3f%10.5f' %
                           (self.zeeman["NOFIELD"][0],
                           self.species,self.expot_lo,
                           self.zeeman["NOFIELD"][1]))
                        if not(self.VdW):
                            out.write('%40s'% (' '))
                        else:
                            out.write('%10.3f%30s' %
                                    (self.VdW, ' '))
                        if not(self.radiative):
                            out.write('%10.3s'% (' '))
                        else:
                            out.write('%10.3f' %
                                    (self.radiative))
                        if not(self.stark):
                            out.write('%10s\n'% (' '))
                        else:
                            out.write('%10.3f\n' %
                                    (self.stark))
                    else:
                        out.write('%10.3f%10.5f%10.3f%10.3f' %
                                (self.wl, self.species, self.expot_lo,self.loggf))
                        if not(self.VdW):
                            out.write('%10s%10.3f%20s' %
                                    (' ',self.DissE, ' '))
                        else:
                            out.write('%10.3f%10.3f%20s' %
                                    (self.VdW, self.DissE, ' '))
                        if not(self.radiative):
                            out.write('%10.3s'% (' '))
                        else:
                            out.write('%10.3f' %
                                    (self.radiative))
                        if not(self.stark):
                            out.write('%10s\n'% (' '))
                        else:
                            out.write('%10.3f\n' %
                                    (self.stark))
    def zeeman_splitting(self, B=None, **kwargs):
        if B:
            self.Bfield = B
        
        self.zeeman["NOFIELD"] = [self.wl, self.loggf]
        self.compute_zeeman_transitions(**kwargs)
        wl = []
        lgf = []
        for transition in self.pi_transitions:
            if (transition.weight > 0):
                wl.append(transition.wavelength)
                lgf.append(numpy.log10(transition.weight*
                    10.0**(self.loggf)))
        self.zeeman["PI"] = [numpy.array(wl), numpy.array(lgf)]

        wl = []
        lgf = []
        for transition in self.lcp_transitions:
            if (transition.weight > 0):
                wl.append(transition.wavelength)
                lgf.append(numpy.log10(transition.weight*
                    10.0**(self.loggf)))
        self.zeeman["LCP"] = [numpy.array(wl), numpy.array(lgf)]

        wl = []
        lgf = []
        for transition in self.rcp_transitions:
            if (transition.weight > 0):
                wl.append(transition.wavelength)
                lgf.append(numpy.log10(transition.weight*
                    10.0**(self.loggf)))
        self.zeeman["RCP"] = [numpy.array(wl), numpy.array(lgf)]

    def compute_zeeman_transitions(self, **kwargs):

        B = self.Bfield

        bohr_magneton = 5.78838176e-5        #eV*T^-1
        hc = 12400                           #eV*Angstroms
        lower_energies = {}
        upper_energies = {}
        for mj in self.lower.mj:
            lower_energies[mj]=self.lower.E+mj*self.lower.g*bohr_magneton*B

        for mj in self.upper.mj:
            upper_energies[mj] = self.upper.E+mj*self.upper.g*bohr_magneton*B

        pi_transitions = []
        lcp_transitions = []
        rcp_transitions = []

        pi_weight = 0.0
        lcp_weight = 0.0
        rcp_weight = 0.0

        delta_J = self.upper.J - self.lower.J
        J1 = self.lower.J

        self.geff = (0.5*(self.lower.g+self.upper.g)
                +0.25*(self.lower.g-self.upper.g)*(self.lower.J*(self.lower.J+1)-
                self.upper.J*(self.upper.J+1.0)))

        for mj in lower_energies.keys():
            if (delta_J == 0.0):
                if upper_energies.has_key(mj+1.0):    #delta Mj = +1 sigma component
                    weight = (J1-mj)*(J1+mj+1.0)
                    rcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj+1]-lower_energies[mj]), weight,
                        mj+1, mj))
                    rcp_weight+=weight
                if upper_energies.has_key(mj):    #delta Mj = 0 Pi component
                    weight = mj**2.0
                    pi_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj]-lower_energies[mj]), weight,
                        mj, mj))
                    pi_weight+=weight
                if upper_energies.has_key(mj-1.0):    #delta Mj = -1 sigma component
                    weight = (J1+mj)*(J1-mj+1.0)
                    lcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj-1]-lower_energies[mj]), weight,
                        mj-1, mj))
                    lcp_weight+=weight
            if (delta_J == 1.0):
                if upper_energies.has_key(mj+1.0):    #delta Mj = +1 sigma component
                    weight = (J1+mj+1.0)*(J1+mj+2.0)
                    rcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj+1]-lower_energies[mj]), weight,
                        mj+1, mj))
                    rcp_weight+=weight
                if upper_energies.has_key(mj):    #delta Mj = 0 Pi component
                    weight = (J1+1.0)**2.0 - mj**2.0
                    pi_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj]-lower_energies[mj]), weight,
                        mj, mj))
                    pi_weight+=weight
                if upper_energies.has_key(mj-1.0):    #delta Mj = -1 sigma component
                    weight = (J1-mj+1.0)*(J1-mj+2.0)
                    lcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj-1]-lower_energies[mj]), weight,
                        mj-1, mj))
                    lcp_weight+=weight
            if (delta_J == -1.0):
                if upper_energies.has_key(mj+1.0):    #delta Mj = +1 sigma component
                    weight = (J1-mj)*(J1-mj-1.0)
                    rcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj+1]-lower_energies[mj]), weight,
                        mj+1, mj))
                    rcp_weight+=weight
                if upper_energies.has_key(mj):    #delta Mj = 0 Pi component
                    weight = J1**2.0 - mj**2.0
                    pi_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj]-lower_energies[mj]), weight,
                        mj, mj))
                    pi_weight+=weight
                if upper_energies.has_key(mj-1.0):    #delta Mj = -1 sigma component
                    weight = (J1+mj)*(J1+mj-1.0)
                    lcp_transitions.append(zeemanTransition(hc/
                        (upper_energies[mj-1]-lower_energies[mj]), weight,
                        mj-1, mj))
                    lcp_weight+=weight

        for transition in rcp_transitions:
            transition.weight /= rcp_weight
        for transition in lcp_transitions:
            transition.weight /= lcp_weight
        for transition in pi_transitions:
            transition.weight /= pi_weight

        self.pi_transitions = pi_transitions
        self.lcp_transitions = lcp_transitions
        self.rcp_transitions = rcp_transitions

class dummy_Line( Spectral_Line ):
    def __init__(self):
        self.wl = 0.0
        self.species = 26.0
        self.element = 26.0
        self.ionization = 0.0
        self.loggf = -5.5
        self.expot_lo = 7.8
        self.Bfield = 0
        self.VdW = None
        self.radiative =None
        self.stark = None
        self.DissE = None
        self.zeeman = {}
        self.zeeman["NOFIELD"] = [self.wl,self.loggf]

    def create(self, wl, outfile):
        self.wl = wl
        self.zeeman["NOFIELD"][0] = wl
        self.dump(out=outfile, mode='MOOGSCALAR')

class MOOG_Line( Spectral_Line ):
    def __init__(self, line, **kwargs):
        self.wl = float(line[0:11])
        self.species = float(line[10:21])
        self.element = numpy.round(self.species)
        self.ionization = (self.species - self.element)*10.0
        self.loggf = float(line[30:41])
        self.loggfHistory = []
        self.VdWHistory = []
        self.expot_lo = float(line[20:31])
        self.Bfield = 0.0
        self.zeeman = {}
        self.zeeman["NOFIELD"] = [self.wl, self.loggf]
        try:
            self.VdW = float(line[40:51])
        except:
            self.VdW = None
        try:
            self.radiative = float(line[80:91])
        except:
            self.radiative = None
        try:
            self.stark = float(line[90:101])
        except:
            self.stark = None
        try:
            self.DissE = float(line[50:61])
        except:
            self.DissE = None

class VALD_Line( Spectral_Line ):
    def __init__(self, line1, line2='', pt='', **kwargs):
        l1 = line1.split(',')
        l2 = line2.split()
        self.element = pt.translate(l1[0].strip('\'').split()[0])
        self.ionization = int(l1[0].strip('\'').split()[1])-1
        self.species = self.element + self.ionization/10.0
        self.wl = float(l1[1])
        self.loggf = float(l1[2])
        self.expot_lo = float(l1[3])
        self.J_lo = float(l1[4])
        self.expot_hi = float(l1[5])
        self.J_hi = float(l1[6])
        self.g_lo = float(l1[7])
        self.g_hi = float(l1[8])
        self.g_eff = float(l1[9])
        self.radiative = float(l1[10])
        self.stark = float(l1[11])
        self.VdW = float(l1[12])
        self.loggfHistory = []
        self.VdWHistory = []
        self.DissE = None
        self.transition = line2.strip().strip('\'')
        self.Bfield = 0.0

        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]
        else:
            self.verbose = False

        if (self.g_lo == 99.0):
            if not (self.species in [70.1, 25.2]):
                angmom = {"S":0, "P":1, "D":2, "F":3, "G":4, "H":5,
                        "I":6, "K":7, "L":8, "M":9}
                n = 0
                try:
                    for char in self.transition:
                        if char.isdigit():
                            S = (float(char)-1.0)/2.0
                        if ((char.isupper()) & (n < 2)):
                            n+=1
                            L = angmom[char]
                            if n == 1:
                                if (self.J_lo > 0.0):
                                    self.g_lo = (1.5+(S*(S+1.0)-L*(L+1))/
                                            (2*self.J_lo*(self.J_lo+1)))
                                else:
                                    self.g_lo = 0.0
                            else:
                                if (self.J_hi > 0.0):
                                    self.g_hi = (1.5+(S*(S+1.0)-L*(L+1))/
                                            (2*self.J_hi*(self.J_hi+1)))
                                else:
                                    self.g_hi = 0.0
                except:
                    self.g_lo = 0.0
                    self.g_hi = 0.0
                    if self.verbose:
                        print("Parsing VALD Transition Failed! %f" % self.wl)
                        print("%s\n" % self.transition)
            else:
                self.g_lo = 0.0
                self.g_hi = 0.0
   
        self.lower = Energy_Level(self.J_lo, self.g_lo, self.expot_lo)
        self.upper = Energy_Level(self.J_hi, self.g_hi,
                self.expot_lo+12400.0/self.wl)
        self.zeeman = {}
        self.zeeman["NOFIELD"] = [self.wl,self.loggf]

class New_VALD_Line( Spectral_Line ):
    def __init__(self, lines, pt='', **kwargs):
        l1 = lines[0].split(',')
        self.ID = l1[0].strip('\'').split()
        self.element = pt.translate(self.ID[0])
        self.Bfield = 0.0

        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]
        else:
            self.verbose = False

        if self.element > 1000:
            isotopes = lines[3].strip().strip('\'').split()[-1].split('(')
            A1 = int(isotopes[1].split(')')[0])
            if len(isotopes) == 2:
                A2 = 1
            else:
                A2 = int(isotopes[2].split(')')[0])
            heavy = max(A1, A2)
            light = min(A1, A2)
            self.species = numpy.float("%4.1f%02d%02d" % (self.element/10., light, heavy))
            #self.species = "%4.1f%2d%2d" % (self.element/10., (10.0+float(l[10])/10)+0.0001*(10.0+float(l[10])%10)
            self.wl = float(l1[1])
            self.loggf = float(l1[2])
            self.expot_lo = float(l1[3])
            self.J_lo = float(l1[4])
            self.expot_hi = float(l1[5])
            self.J_hi = float(l1[6])
            self.g_lo = float(l1[7])
            self.g_hi = float(l1[8])
            self.g_eff = float(l1[9])
            self.radiative = float(l1[10])
            self.stark = float(l1[11])
            self.VdW = float(l1[12])
            self.loggfHistory = []
            self.VdWHistory = []
            self.DissE = pt.DissE(self.ID[0])
        else:
            self.ionization = int(self.ID[1])-1
            self.species = self.element + self.ionization/10.0
            self.wl = float(l1[1])             # WL in Air
            self.loggf = float(l1[2])
            self.expot_lo = float(l1[3])       # in eV
            self.J_lo = float(l1[4])
            self.expot_hi = float(l1[5])
            self.J_hi = float(l1[6])
            self.g_lo = float(l1[7])
            self.g_hi = float(l1[8])
            self.g_eff = float(l1[9])
            self.radiative = float(l1[10])
            if self.element == 1.0:
                self.stark = -4.0
                self.VdW = -5.4
            else:
                self.stark = float(l1[11])
                self.VdW = float(l1[12])
            self.loggfHistory = []
            self.VdWHistory = []
            self.DissE = None
            self.upperTerm = lines[1].strip().strip('\'').split()
            self.lowerTerm = lines[2].strip().strip('\'').split()
            self.references = lines[3].strip().strip('\'')

            if (self.g_lo == 99.0):
                # Lower State
                try:
                    if self.lowerTerm[0] == 'LS':
                        self.g_lo = self.parse_LS_coupling(self.lowerTerm[-1], self.J_lo)
                    elif self.lowerTerm[0] == 'JJ':
                        #print("Awww Shucks, the JJ coupling parser isn't ready yet!")
                        #self.g_lo = self.parse_JJ_coupling(self.lowerTerm[-1])
                        self.g_lo = 0.0
                    elif self.lowerTerm[0] == 'JK':
                        #print("Awww Shucks, the JK coupling parser isn't ready yet!")
                        self.g_lo = 0.0
                    elif self.lowerTerm[0] == 'LK':
                        #print("Awww Shucks, the LK coupling parser isn't ready yet!")
                        self.g_lo = 0.0
                    else:
                        self.g_lo = 0.0
                except:
                    print("Lower state of line at %.3f failed!" % self.wl)
                    self.g_lo = 0.0
    
                # Upper State
                try:
                    if self.upperTerm[0] == 'LS':
                        self.g_hi = self.parse_LS_coupling(self.upperTerm[-1], self.J_hi)
                    elif self.upperTerm[0] == 'JJ':
                        #print("Awww Shucks, the JJ coupling parser isn't ready yet!")
                        self.g_hi = 0.0
                    elif self.upperTerm[0] == 'JK':
                        #print("Awww Shucks, the JK coupling parser isn't ready yet!")
                        self.g_hi = 0.0
                    elif self.upperTerm[0] == 'LK':
                        #print("Awww Shucks, the LK coupling parser isn't ready yet!")
                        self.g_hi = 0.0
                    else:
                        self.g_hi = 0.0
                except:
                    print("Upper state of line at %.3f failed!" % self.wl)
                    self.g_hi = 0.0

        self.lower = Energy_Level(self.J_lo, self.g_lo, self.expot_lo)
        self.upper = Energy_Level(self.J_hi, self.g_hi,
                self.expot_lo+12400.0/self.wl)
        self.zeeman = {}
        self.zeeman["NOFIELD"] = [self.wl,self.loggf]

    def parse_LS_coupling(self, term, J):
        if J <= 0:
            return 0.0
        angmom = {"S":0, "P":1, "D":2, "F":3, "G":4, "H":5,
                "I":6, "K":7, "L":8, "M":9}
        try:
            if term[0].isdigit():
                S = (float(term[0])-1.0)/2.0
            else:
                S = (float(term[1])-1.0)/2.0
            if term[-1] == '*':
                L = angmom[term[-2]]
            else:
                L = angmom[term[-1]]
            lande_g = (1.5+(S*(S+1.0)-L*(L+1))/
                    (2.0*J*(J+1.0)))
        except:
            print("Parsing of LS term at %.3f FAILED!!!" % self.wl)
            lande_g = 0.0
        return lande_g

    def parse_JJ_coupling(self, term):
        term = term.strip('*').replace('(','').replace(')','').split(',')
        J1 = term[0].split('/')
        if len(J1) == 2:
            J1 = float(J1[0])/float(J1[1])
        else:
            J1 = float(J1[0])
        J2 = term[1].split('/')
        if len(J2) == 2:
            J2 = float(J2[0])/float(J2[1])
        else:
            J2 = float(J2[0])
        lande_g = 0.0
        return lande_g


class Plez_CN_Line( Spectral_Line ):
    def __init__(self, line, **kwargs):
        l = line.split()
        self.wl = float(l[0])
        self.species = float(l[1])
        self.DissE = 7.72
        self.expot_lo = float(l[2])
        self.loggf = float(l[3])
        self.VdW = None
        self.radiative = None
        self.stark = None
        self.zeeman = {}
        self.transition = None
        self.J_lo = None
        self.J_hi = None
        self.g_lo = None
        self.g_hi = None
        self.g_eff = None
        self.loggfHistory = []
        self.VdWHistory = []
        self.zeeman["NOFIELD"] = [self.wl, self.loggf]
        self.Bfield = 0.0

        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]
        else:
            self.verbose = False

class Goorvitch_CO_Line( Spectral_Line ):
    def __init__(self, line, **kwargs):
        l = line.split('|')
        try:
            self.wl = 1.0e8/float(l[0])
            self.species = 0608.0+0.001*(10.0+float(l[10])/10)+0.0001*(10.0+float(l[10])%10)
            self.DissE = 11.10
            self.expot_lo = 1.23981e-4*float(l[3])
            self.loggf = numpy.log10(float(l[4]))
            self.VdW = None
            self.radiative = None
            self.stark = None
            self.zeeman = {}
            self.transition = None
            self.J_lo = None
            self.J_hi = None
            self.g_lo = None
            self.g_hi = None
            self.g_eff = None
            self.loggfHistory = []
            self.VdWHistory = []
            self.zeeman["NOFIELD"] = [self.wl, self.loggf]
            self.Bfield = 0.0

            if "verbose" in kwargs:
                self.verbose = kwargs["verbose"]
            else:
                self.verbose = False

        except:
            self.wl = -99.9

class HITRAN_Line( Spectral_Line ):
    def __init__(self, line, hitran_dictionary, **kwargs):
        hitran_code = int(line[0:2])
        isotope_code = int(line[2])
        self.species = hitran_dictionary.isotopes[hitran_code][isotope_code]
        self.DissE = hitran_dictionary.DissE[hitran_code]
        self.wl = 10000.0/float(line[3:15])*10000.0/1.000273
        self.expot_lo = 1.23986e-4*float(line[45:56])
        Einstein_A = float(line[26:35])
        g_up = float(line[145:154])
        g_low = float(line[154:])
        self.loggf = numpy.log10(1.4991e-16*self.wl**2*g_up*Einstein_A)
        self.VdW = None
        self.radiative = None
        self.stark = None
        self.zeeman = {}
        self.transition = None
        self.J_lo = None
        self.J_hi = None
        self.g_lo = None
        self.g_hi = None
        self.g_eff = None
        self.loggfHistory = []
        self.VdWHistory = []
        self.zeeman["NOFIELD"] = [self.wl, self.loggf]
        self.Bfield = 0.0

        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]
        else:
            self.verbose = False


class zeemanTransition( object):
    def __init__(self, wavelength, weight, m_up, m_low):
        self.wavelength = wavelength
        self.weight = weight
        self.m_up = m_up
        self.m_low = m_low

    def __eq__(self, other):
        return ( (self.wavelength == other.wavelength) &
                (self.m_up == other.m_up) & (self.m_low == other.m_low) )

class Energy_Level( object ):
    def __init__(self, J, g, E):
        self.E = E
        self.J = J
        if g != 99:
            self.g = g
        else:
            self.g = 1.0

        self.mj = numpy.arange(self.J, (-1.0*self.J)-0.5, step = -1)

def generate_CorrectedLines(original_files, new_files, outfile, compfile):
    out = open(outfile, 'w')
    comp = open(compfile, 'w')

    outlines = []
    complines = []

    for orig, new in zip(original_files, new_files):
        with open(orig) as o, open(new) as n:
            old_lines = o.readlines()
            new_lines = n.readlines()

            for ol, nl in zip(old_lines, new_lines):
                if ol != nl:
                    outlines.append(nl)
                    complines.append(ol)

    order = numpy.argsort(outlines)
    out.writelines(numpy.array(outlines)[order])
    comp.writelines(numpy.array(complines)[order])
    out.close()
    comp.close()

def parse_VALD(VALD_list, strong_file, wl_start, wl_stop, Bfield,
        gf_corrections):
    pt = periodicTable()

    corrected = []
    for line in open(gf_corrections, 'r'):
        corrected.append(MOOG_Line(line))

    strong = []
    for line in open(strong_file, 'r'):
        l = line.split()
        strong.append([float(l[0]), float(l[1])])
    
    vald_in = open(VALD_list, 'r')
    l1 = ''
    l2 = ''
    stronglines = []
    weaklines = []
    for line in vald_in:
        if line[0] != '#':
            if line[0] == '\'':
                l1 = line
            else:
                l2 = line
                current_line = VALD_Line(l1, l2, pt)
                wl = current_line.wl
                if ( (wl > wl_start) & (wl < wl_stop) ):
                    for cl in corrected:
                        if (cl.wl == current_line.wl) & (cl.expot_lo ==
                                current_line.expot_lo):
                            print("Making a correction!")
                            current_line.loggf = cl.loggf
                            current_line.zeeman["NOFIELD"][1] = cl.loggf
                            current_line.VdW = cl.VdW
                            current_line.stark = cl.stark
                            current_line.radiative = cl.radiative
                    current_line.zeeman_splitting(Bfield)
                    species = current_line.species
                    if ( [wl, species] in strong):
                        stronglines.append(current_line)
                    else:
                        weaklines.append(current_line)

    return stronglines, weaklines

def parse_HITRAN(HITRAN_file, wl_start, wl_stop, B_field, 
        gf_corrections, **kwargs):

    corrected = []
    for line in open(gf_corrections, 'r'):
        corrected.append(MOOG_Line(line))

    ht = HITRAN_Dictionary()
    hitran_in = open(HITRAN_file, 'r')
    lines = []
    for line in hitran_in:
        current_line = HITRAN_Line(line, ht)
        if ( (current_line.wl > wl_start) & (current_line.wl < wl_stop) ):
            for cl in corrected:
                if (cl.wl == current_line.wl) & (cl.expot_lo ==
                         current_line.expot_lo):
                    current_line.loggf = cl.loggf
                    current_line.VdW = cl.VdW
                    current_line.stark = cl.stark
                    current_line.radiative = cl.radiative
            if "weedout" in kwargs:
                if current_line.expot_lo < kwargs["weedout"]:
                    lines.append(current_line)
                else:
                    print('Tossed CO line!')
            else:
                lines.append(current_line)

    return lines

def parse_Plez_CN(CN_file, wl_start, wl_stop, B_field, gf_corrections, 
        **kwargs):

    corrected = []
    for line in open(gf_corrections, 'r'):
        corrected.append(MOOG_Line(line))

    cn_in = open(CN_file, 'r')
    lines = []
    for line in cn_in:
        current_line = Plez_CN_Line(line)
        if ( (current_line.wl > wl_start) & (current_line.wl < wl_stop) ):
            for cl in corrected:
                if (cl.wl == current_line.wl) & (cl.expot_lo ==
                         current_line.expot_lo):
                    current_line.loggf = cl.loggf
                    current_line.VdW = cl.VdW
                    current_line.stark = cl.stark
                    current_line.radiative = cl.radiative
            lines.append(current_line)

    return lines

def parse_Goorvitch_CO(CO_file, wl_start, wl_stop, B_field, gf_corrections,
        **kwargs):

    corrected = []
    for line in open(gf_corrections, 'r'):
        corrected.append(MOOG_Line(line))

    co_in = open(CO_file, 'r')
    lines = []
    for line in co_in:
        current_line = Goorvitch_CO_Line(line)
        if ( (current_line.wl > wl_start) & (current_line.wl < wl_stop) ):
            for cl in corrected:
                if (cl.wl == current_line.wl) & (cl.expot_lo ==
                         current_line.expot_lo):
                    current_line.loggf = cl.loggf
                    current_line.VdW = cl.VdW
                    current_line.stark = cl.stark
                    current_line.radiative = cl.radiative
            lines.append(current_line)

    return lines


def write_par_file(wl_start, wl_stop, stage_dir, b_dir, prefix, temps=None, 
        gravs=None, mode='gridstokes', strongLines=False, **kwargs):
    if mode=='gridstokes':
        fn = 'batch.par'
        suffix = '.stokes'
    elif mode == 'gridsyn':
        fn = 'batch.gridsyn'
        suffix = '.scalar'
    elif mode == 'stokes':
        fn = 'batch.stokes'
        suffix = '.stokes'
    elif mode == 'synth':
        fn = 'batch.synth'
        suffix = '.scalar'

    outfile_name = os.path.join(stage_dir,'Parfiles', b_dir, fn)

    if "OUT_DIR" in kwargs.keys():
        output_prefix = kwargs["OUT_DIR"]
    else:
        output_prefix = '../../Output/'+b_dir+'/'

    line_prefix = '../../Linelists/'+b_dir+'/'

    labels = {'terminal':'x11',
            'strong':1, 
            'atmosphere':1, 
            'molecules':2,
            'lines':1,
            'damping':1,
            'freeform':2,
            'flux/int':0,
            'diskflag':1}
            #'plot':2, 
            #'obspectrum':5}
    file_labels = {'summary_out':'../../Output/'+b_dir+'/summary.out',
            'standard_out':output_prefix+'out1',
            'smoothed_out':output_prefix+'smoothed.out',
            'atmos_dir':'/home/deen/Data/Atmospheres/MARCS/',
            'out_dir':output_prefix,
            'lines_in':line_prefix+prefix+'_weak_linelist'+suffix,
            'stronglines_in':line_prefix+prefix+'_strong_linelist'+suffix}
            #'model_in':'model.md',
            #'observed_in':'observed.dat'}


    for l in labels:
        if l in kwargs:
            labels[l] = kwargs[l]
            
    for fl in file_labels:
        if fl in kwargs:
            file_labels[fl] = kwargs[fl]

    pf = open(outfile_name, 'w')

    pf.write(mode+'\n')
    for fl in file_labels:
        pf.write(fl+'   \''+file_labels[fl]+'\'\n')
    for l in labels:
        pf.write(l+'    '+str(labels[l])+'\n')

    pf.write('synlimits\n')
    pf.write('               '+str(wl_start)+' '
             +str(wl_stop)+' 0.01 3.50\n')
    pf.write('plotpars       1\n')
    pf.write('               '+str(wl_start)+' '
             +str(wl_stop)+' 0.02 1.00\n')
    pf.write('               0.00 0.000 0.000 1.00\n')
    pf.write('               g 0.150 0.00 0.00 0.00 0.00\n')

    if ( (mode=='gridstokes') | (mode=='gridsyn')):
        run_number = 1

        if (not temps):
            temps = range(2500, 4100, 100)+range(4250, 6250, 250)
        if (not gravs):
            gravs = range(300, 550, 50)

        for T in temps:
            for G in gravs:
                pf.write('RUN            '+str(run_number)+'\n')
                if mode == 'gridstokes':
                    pf.write('stokes_out   \''+prefix+
                        '_MARCS_T'+str(T)+'G'+str(G)+'\'\n')
                else:
                    pf.write('smoothed_out   \''+prefix+
                        '_MARCS_T'+str(T)+'G'+str(G)+'\'\n')
                pf.write('hardpost_out   \'../../Output/'+b_dir+'/dummy.ps\'\n')
                pf.write('model_in       \'MARCS_T'+
                        str(T)+'_G'+str(G/100.0)+'_M0.0_t2.0.md\'\n')
                pf.write('abundances     1  1\n')
                pf.write('    12      0.0\n')
                run_number += 1
    pf.close()


class Angle( object ):
    def __init__(self, line):
        l = line.split()
        self.n = int(l[0])
        self.az = [float(l[1]), float(l[2]), float(l[3])]
        self.longitude = [float(l[4]), float(l[5])]
        self.phi = float(l[6])
        self.chi = float(l[7])
        self.mu = float(l[8])

class Diskoball( object ):
    def __init__(self, name, **kwargs):
        self.name = name
        if "DIR" in kwargs.keys():
            self.directory = kwargs["DIR"]
        else:
            self.directory = '../'
        if "VSINI" in kwargs.keys():
            self.vsini = kwargs["VSINI"]
        else:
            self.vsini = 0.0
        self.dfI = self.directory+self.name+'.spectrum_I'
        self.dfQ = self.directory+self.name+'.spectrum_Q'
        self.dfU = self.directory+self.name+'.spectrum_U'
        self.dfV = self.directory+self.name+'.spectrum_V'
        self.dfCont = self.directory+self.name+'.continuum'
        self.dfAngles = self.directory+self.name+'.angles'

        Angles = open(self.dfAngles, 'r')
        StokesI = open(self.dfI, 'r')
        StokesQ = open(self.dfQ, 'r')
        StokesU = open(self.dfU, 'r')
        StokesV = open(self.dfV, 'r')
        Continuum = open(self.dfCont, 'r')

        linecounter = 0
        self.ang_info = []
        for line in Angles:
            if linecounter == 0:
                l = line.split()
                self.ncells = int(l[0])
                self.nrings = int(l[1])
                self.inclination = float(l[2])
                self.PA = float(l[3])
                self.cell_area = 4.0*3.1415926/self.ncells
                linecounter +=1
            else:
                self.ang_info.append(Angle(line))


        wl = []
        I = []
        Q = []
        U = []
        V = []
        C = []
        
        for line in StokesI:
            l = line.split()
            wl.append(float(l[0]))
            a = []
            for fluxes in l[1:]:
                try:
                    a.append(float(fluxes))
                except:
                    print("Warning! I crazy format")
                    a.append(float(0.0))
            
            I.append(a)

        for line in StokesQ:
            l = line.split()
            a = []
            for fluxes in l[1:]:
                try:
                    a.append(float(fluxes))
                except:
                    print("Warning! Q crazy format")
                    a.append(float(0.0))

            Q.append(a)

        for line in StokesU:
            l = line.split()
            a = []
            for fluxes in l[1:]:
                try:
                    a.append(float(fluxes))
                except:
                    print("Warning! U crazy format")
                    a.append(float(0.0))

            U.append(a)

        for line in StokesV:
            l = line.split()
            a = []
            for fluxes in l[1:]:
                try:
                    a.append(float(fluxes))
                except:
                    print("Warning! V crazy format")
                    a.append(float(0.0))

            V.append(a)

        for line in Continuum:
            l = line.split()
            a = []
            for fluxes in l[1:]:
                try:
                    a.append(float(fluxes))
                except:
                    print("Warning! C crazy format")
                    a.append(float(0.0))

            C.append(a)
        
        self.wl = numpy.array(wl)
        I = numpy.array(I)
        Q = numpy.array(Q)
        U = numpy.array(U)
        V = numpy.array(V)
        C = numpy.array(C)
        self.I = I.transpose()
        self.Q = Q.transpose()
        self.U = U.transpose()
        self.V = V.transpose()
        self.C = C.transpose()

        wave = numpy.mean(self.wl)
        if ((1.0/(wave/10000.0)) < 2.4):
            self.alpha = -0.023 + 0.292/(wave/10000.0)
        else:
            self.alpha = -0.507 + 0.441/(wave/10000.0)


    def interpolate(self, stepsize):
        self.wave = numpy.arange(self.wl[0], self.wl[-1], step=stepsize)
        fI = scipy.interpolate.UnivariateSpline(self.wl, self.integrated_I, s=0)
        fQ = scipy.interpolate.UnivariateSpline(self.wl, self.integrated_Q, s=0)
        fU = scipy.interpolate.UnivariateSpline(self.wl, self.integrated_U, s=0)
        fV = scipy.interpolate.UnivariateSpline(self.wl, self.integrated_V, s=0)
        fC = scipy.interpolate.UnivariateSpline(self.wl, self.integrated_C, s=0)
        self.flux_I = fI(self.wave)
        self.flux_Q = fQ(self.wave)
        self.flux_U = fU(self.wave)
        self.flux_V = fV(self.wave)
        self.flux_C = fC(self.wave)


    def disko(self):
        r2d = 180.0/numpy.pi
        final_I = numpy.zeros(len(self.wl))
        final_Q = numpy.zeros(len(self.wl))
        final_U = numpy.zeros(len(self.wl))
        final_V = numpy.zeros(len(self.wl))
        final_C = numpy.zeros(len(self.wl))
        
        total_weight = 0.0
        T_I = numpy.matrix([[1.0, 0.0, 0.0],
              [0.0, numpy.cos(self.inclination), numpy.sin(self.inclination)],
              [0.0, -numpy.sin(self.inclination), numpy.cos(self.inclination)]])

        emergent_vector = numpy.matrix([1.0, 0.0, 0.0])

        for tile in zip(self.I, self.Q, self.U, self.V, self.C, self.ang_info):
            azimuth = tile[5].az
            n_az_steps = int(azimuth[2]*r2d-azimuth[1]*r2d)
            azs = azimuth[1]+(numpy.arange(n_az_steps)+0.5)*(azimuth[2]-
                    azimuth[1])/n_az_steps
            az1 = azimuth[1]+(numpy.arange(n_az_steps))*(azimuth[2]-
                    azimuth[1])/n_az_steps
            az2 = azimuth[1]+(numpy.arange(n_az_steps)+1.0)*(azimuth[2]-
                    azimuth[1])/n_az_steps
            longitude = tile[5].longitude
            dphi = longitude[1]
            n_phi_steps = int(dphi*r2d)
            phis = longitude[0]-dphi/2.0+(numpy.arange(n_phi_steps)+
                    0.5)*dphi/n_phi_steps
            for az in zip(azs, az1, az2):
                T_rho = numpy.matrix([[0.0, 0.0, 1.0],
                        [-numpy.cos(az[0]), numpy.sin(az[0]), 0.0],
                        [numpy.sin(az[0]), numpy.cos(az[0]), 0.0]])
                daz = numpy.sin(az[2])-numpy.sin(az[1])
                area = daz*dphi/n_phi_steps
                for phi in phis:
                    T_eta = numpy.matrix([
                            [numpy.cos(phi), -numpy.sin(phi), 0.0],
                            [numpy.sin(phi), numpy.cos(phi), 0.0],
                            [0.0, 0.0, 1.0]])
                    surface_vector = T_I*T_eta*T_rho*emergent_vector.T
                    mu = surface_vector.A[2][0]
                    if (mu > 0.00001):
                        projected_area = area*mu#/(4.0*pi)
                        limb_darkening = (1.0-(1.0-mu**self.alpha))
                        weight = projected_area*limb_darkening
                        total_weight += weight
                        final_I = final_I + weight*tile[0]/tile[4]
                        final_Q = final_Q + weight*tile[1]/tile[4]
                        final_U = final_U + weight*tile[2]/tile[4]
                        final_V = final_V + weight*tile[3]/tile[4]
                        final_C = final_C + weight*tile[4]

        self.integrated_I = final_I/total_weight
        self.integrated_Q = final_Q/total_weight
        self.integrated_U = final_U/total_weight
        self.integrated_V = final_V/total_weight
        self.integrated_C = final_C/total_weight

    def save(self, outfile):
        SpectralTools.write_2col_spectrum(outfile+'.I', self.wave, self.flux_I)
        SpectralTools.write_2col_spectrum(outfile+'.Q', self.wave, self.flux_Q)
        SpectralTools.write_2col_spectrum(outfile+'.U', self.wave, self.flux_U)
        SpectralTools.write_2col_spectrum(outfile+'.V', self.wave, self.flux_V)
        SpectralTools.write_2col_spectrum(outfile+'.C', self.wave, self.flux_C)

class MoogStokes_IV_Spectrum( object ):
    #"""
    def __init__(self, name='', memory=False, **kwargs):
        self.memory = memory
        if self.memory:
            self.parent = kwargs["PARENT"]
            self.deltav = self.parent.deltav
            self.vsini = self.parent.vsini
        else:
            self.name = name
            if "DIR" in kwargs.keys():
                self.directory = kwargs["DIR"]
            else:
                self.directory = '../'
            if "DELTAV" in kwargs.keys():
                self.deltav = kwargs["DELTAV"]
            else:
                self.deltav = 0.1            #  wl spacing in km/s
            if "VSINI" in kwargs.keys():
                self.vsini = kwargs["VSINI"]
            else:
                self.vsini = 0.0

            self.angle_file = self.directory+self.name+'.angles'
            self.continuum_file = self.directory+self.name+'.continuum'
            self.I_file = self.directory+self.name+'.spectrum_I'
            self.V_file = self.directory+self.name+'.spectrum_V'

            self.nangles = 0
            self.phi = []
            self.mu = []
            self.wl = []
            self.I = []
            self.V = []
            self.continuum = []

        self.loadAngles()
        self.loadSpectra()

        self.interpolateSpectra()
        self.diskInt()
    #"""


    def loadAngles(self):
        if self.memory:
            self.phi = self.parent.phi_angle[:self.parent.ncells]
            self.mu = self.parent.mus[:self.parent.ncells]
        else:
            df = open(self.angle_file, 'r')
            for line in df:
                l = line.split()
                if len(l) == 1:
                    self.nangles = int(l[0])
                else:
                    self.phi.append(float(l[1]))
                    self.mu.append(float(l[2]))

            self.phi = numpy.array(self.phi)
            self.mu = numpy.array(self.mu)

    def loadSpectra(self):
        if self.memory:
            self.I = numpy.array(self.parent.flux_I)/numpy.array(self.parent.continuum)
            self.V = numpy.array(self.parent.flux_V)/numpy.array(self.parent.continuum)
            self.continuum = numpy.array(self.parent.continuum)
            self.wl = numpy.array(self.parent.wave)
        else:
            df_I = open(self.I_file, 'r')
            df_V = open(self.V_file, 'r')
            df_C = open(self.continuum_file, 'r')

            continuum = []
            I = []
            V = []
            wl = []
            for line in df_C:
                l = line.split()
                wl.append(float(l[0]))
                a = []
                for fluxes in l[1:]:
                    try:
                        a.append(float(fluxes))
                    except:
                        print("Warning! Crazy Continuum format!", fluxes)
                        a.append(float(0.0))
                continuum.append(a)

            for line in df_I:
                l = line.split()
                a = []
                for fluxes in l[1:]:
                    try:
                        a.append(float(fluxes))
                    except:
                        print("Woah there pardner! Crazy format - Stokes I!", fluxes)
                        a.append(float(0.0))
                I.append(a)

            for line in df_V:
                l = line.split()
                a = []
                for fluxes in l[1:]:
                    try:
                        a.append(float(fluxes))
                    except:
                        print("Woah there pardner! Crazy format - Stokes V!", fluxes)
                        a.append(float(0.0))
                V.append(a)

            self.wl = numpy.array(wl)
            I = numpy.array(I)
            V = numpy.array(V)
            continuum = numpy.array(continuum)
            self.continuum = continuum.transpose()
            self.I = I.transpose()/self.continuum
            self.V = V.transpose()/self.continuum

        wave = numpy.mean(self.wl)
        if ((1.0/(wave/10000.0)) < 2.4):
            self.alpha = -0.023 + 0.292/(wave/10000.0)
        else:
            self.alpha = -0.507 + 0.441/(wave/10000.0)

    def interpolateSpectra(self):
        deltav = self.deltav
        c = 3e5                        #km/s
        wl_start = numpy.min(self.wl)
        wl_max = numpy.max(self.wl)
        new_wl = []
        new_wl.append(wl_start)
        while new_wl[-1] < wl_max:
            d_lambda = new_wl[-1]*deltav/c
            new_wl.append(new_wl[-1]+d_lambda)
        self.new_wl = numpy.array(new_wl[0:-1])

        new_I = []
        new_V = []
        for I,V in zip(self.I, self.V):
            fI = scipy.interpolate.UnivariateSpline(self.wl, I, s=0)
            fV = scipy.interpolate.UnivariateSpline(self.wl, V, s=0)
            new_I.append(fI(self.new_wl))
            new_V.append(fV(self.new_wl))

        self.new_I = numpy.array(new_I)
        self.new_V = numpy.array(new_V)

    def diskInt(self):
        deltav = self.deltav
        vsini = self.vsini
        c = 3e5
        limb_darkening = []
        for i in range(len(self.mu)):
            limb_darkening.append(1.0-(1.0-self.mu[i]**(self.alpha)))

        self.limb_darkening = numpy.array(limb_darkening)
        continuum = []
        for i in range(len(self.mu)):
            self.new_I[i] *= self.limb_darkening[i]
            self.new_V[i] *= self.limb_darkening[i]
            continuum.append(numpy.ones(len(self.new_I[i]))
                    *self.limb_darkening[i])

        continuum = numpy.array(continuum)

        self.final_spectrum_I = self.rtint(self.mu, self.new_I,
                continuum, deltav, vsini, 0.0)
        self.final_spectrum_V = self.rtint(self.mu, self.new_V,
                continuum, deltav, vsini, 0.0)
        
    def save(self, outfile):
        SpectralTools.write_2col_spectrum(outfile, self.new_wl, self.final_spectrum_I)

    def rtint(self, mu, inten, cont, deltav, vsini_in, vrt_in, **kwargs):
        """
    This is a python translation of Jeff Valenti's disk integration routine
    
    PURPOSE:
        Produces a flux profile by integrating intensity profiles (sampled
           at various mu angles) over the visible stellar surface.

    Calling Sequence:
        flux = rtint(mu, inten, deltav, vsini, vrt)

    INPUTS:
        MU: list of length nmu cosine of the angle between the outward normal
            and the line of sight for each intensity spectrum INTEN
        INTEN:  list (of length nmu) numpy arrays (each of length npts)
            intensity spectra at specified values of MU
        DELTAV: (scalar) velocity spacing between adjacent spectrum points in
            INTEN (same units as VSINI and VRT)

        VSIN (scalar) maximum radial velocity, due to solid-body rotation
        VRT (scalar) radial-tangential macroturbulence parameter, i.e.. sqrt(2)
            times the standard deviation of a Gaussian distribution of 
            turbulent velocities.  The same distribution function describes
            the raidal motions of one component and the tangential motions of
            a second component.  Each component covers half the stellar surface.
            See "Observation and Analysis of Stellar Photospheres" by Gray.

    INPUT KEYWORDS:
        OSAMP: (scalar) internal oversamping factor for the convolutions.  By
            default, convolutions are done using the input points (OSAMP=1), 
            but when OSAMP is set to higher integer values, the input spectra
            are first oversamping via cubic spline interpolation.

    OUTPUTS:
        function value: numpy array of length npts producing the disk-integrated
            flux profile.

    RESTRICTIONS:
        Intensity profiles are weighted by the fraction of the projected stellar
            surface they represent, apportioning the area between adjacent MU
            points equally.  Additional weights (such as those used in a Gauss-
            Legendre quadrature) cannot meaningfully be used in this scheme.
            About twice as many points are required with this scheme to achieve
            the same precision of Gauss-Legendre quadrature.
        DELTAV, VSINI, and VRT must all be in the same units (e.q. km/s).
        If specified, OSAMP should be a positive integer

    AUTHOR'S REQUEST:
        If you use this algorithm in work that you publish, please cite...

    MODIFICATION HISTORY:
            Feb 88  GM Created ANA version
         13 Oct 92 JAV Adapted from G. Marcy's ANA routine of same name
         03 Nov 93 JAV Switched to annular convolution technique
         12 Nov 93 JAV Fixed bug. Intensity components not added when vsini=0
         14 Jun 94 JAV Reformatted for "public" release.  Heavily commented.
                 Pass deltav instead of 2.998d5/deltav.  Added osamp
                    keyword.  Added rebinning logic and end of routine.
                 Changed default osamp from 3 to 1.
         20 Feb 95 JAV Added mu as an argument to handle arbitrary mu sampling
                    and remove ambiguity in intensity profile ordering.
                 Interpret VTURB as sqrt(2)*sigma instead of just sigma
                 Replaced call_external with call to spl_{init|interp}.
         03 Apr 95 JAV Multiply flux by !pi to give observed flux.
         24 Oct 95 JAV Force "nmk" padding to be at least 3 pixels
         18 Dec 95 JAV Renamed from dkint() to rtint().  No longer make local
                    copy of intensities.  Use radial-tangential instead of 
                    isotropic Gaussian macroturbulence.
         26 Jan 99 JAV For NMU=1 and VSINI=0, assume resolved solar surface;
                    apply R-T macro, but supress vsini broadening.
         01 Apr 99 GMH Use annuli weights, rather than assuming equal area.
         27 Feb 13 CPD Translated to Python

        """
    
        #make local copies of various input vars, which will be altered below
        vsini = float(vsini_in)
        vrt = float(vrt_in)

        if "OSAMP" in kwargs:
            os = max(round(kwargs["OSAMP"]), 1)
        else:
            os = 1

        #Convert input MU to proj. radii, R of annuli for star of unit radius
        #(which is just sine rather than cosine of the angle between the outward
        #normal and the LOS)
        rmu = numpy.sqrt(1.0-mu**2)

        #Sort the proj. radii and corresponding intensity spectra into ascending
        #order (i.e. from disk center to limb), which is equivalent to sorting
        #MU in decending order
        order = numpy.argsort(rmu)
        rmu = rmu[order]
        nmu = len(mu)
        if (nmu == 1):
            vsini = 0.0

        #Calculate the proj. radii for boundaries of disk integration annuli.
        #The n+1 boundaries are selected so that r(i+1) exactly bisects the area
        #between rmu(i) and rmu(i+1).  The innermost boundary, r(0) is set to 0
        #(Disk center) and the outermost boundary r(nmu) is set to to 1 (limb).
        if ((nmu > 1) | (vsini != 0)):
            r = numpy.sqrt(0.5*(rmu[0:-1]**2.0+rmu[1:]**2.0)).tolist()
            r.insert(0, 0.0)
            r.append(1.0)
            r = numpy.array(r)
    
        #Calculate integration weights for each disk integration annulus.  The
        #weight is given by the relative area of each annulus, normalized such
        #that the sum of all weights is unity.  Weights for limb darkening are
        #included explicitly in intensity profiles, so they aren't needed here.
            wt = r[1:]**2.0 - r[0:-1]**2.0
        else:
            wt = numpy.array([1.0])
        
        #Generate index vectors for input and oversampled points.  Note that the
        #oversampled indicies are carefully chosen such that every "os" finely
        #sampled points fit exactly into one input bin.  This makes it simple to
        #"integrate" the finely sampled points at the end of the routine.

        npts = len(inten[0])
        xpix = numpy.arange(npts)
        nfine = os*npts
        xfine = 0.5/os * 2.0*numpy.arange(nfine)-os+1

        #Loop through annuli, constructing and convolving with rotation kernels.
        dummy = 0
        yfine = numpy.zeros(nfine)
        cfine = numpy.zeros(nfine)
        flux = numpy.zeros(nfine)
        continuum = numpy.zeros(nfine)
        for m, y, c, w, i in zip(mu, inten, cont, wt, range(nmu)):
            #use cubic spline routine to make an oversampled version of the
            #intensity profile for the current annulus.
            if os== 1:
                yfine = y.copy()
                cfine = c.copy()
            else:
                yspl = scipy.interpolate.splrep(xpix, y)
                cspl = scipy.interpolate.splref(xpix, c)
                yfine = scipy.interpolate.splev(yspl, xfine)
                cfine = scipy.interpolate.splev(cspl, xfine)

        # Construct the convolution kernel which describes the distribution of 
        # rotational velocities present in the current annulus. The distribution
        # has been derived analyitically for annuli of arbitrary thickness in a 
        # rigidly rotating star.  The kernel is constructed in two places: one 
        # piece for radial velocities less than the maximum velocity along the
        # inner edge of annulus, and one piece for velocities greater than this
        # limit.
            if vsini > 0:
                r1 = r[i]
                r2 = r[i+1]
                dv = deltav/os
                maxv = vsini * r2
                nrk = 2*long(maxv/dv) + 3
                v = dv * (numpy.arange(nrk) - ((nrk-1)/2.))
                rkern = numpy.zeros(nrk)
                j1 = scipy.where(abs(v) < vsini*r1)
                if len(j1[0]) > 0:
                    rkern[j1] = (numpy.sqrt((vsini*r2)**2 - v[j1]**2)-
                            numpy.sqrt((vsini*r1)**2 - v[j1]**2))
                j2 = scipy.where((abs(v) >= vsini*r1) & (abs(v) <= vsini*r2))
                if len(j2[0]) > 0:
                    rkern[j2] = numpy.sqrt((vsini*r2)**2 - v[j2]**2)
                rkern = rkern / rkern.sum()   # normalize kernel


        # Convolve the intensity profile with the rotational velocity kernel for
        # this annulus.  Pad end of each profile with as many points as are in
        # the convolution kernel, reducing Fourier ringing.  The convolution 
        # may also be done with a routine called "externally" which efficiently
        # shifts and adds.
                if nrk > 3:
                    yfine = scipy.convolve(yfine, rkern, mode='same')
                    cfine = scipy.convolve(cfine, rkern, mode='same')

        # Calc projected simga for radial and tangential velocity distributions.
            sigma = os*vrt/numpy.sqrt(2.0) /deltav
            sigr = sigma * m
            sigt = sigma * numpy.sqrt(1.0 - m**2.)

        # Figure out how many points to use in macroturbulence kernel
            nmk = max(min(round(sigma*10), (nfine-3)/2), 3)

        # Construct radial macroturbulence kernel w/ sigma of mu*VRT/sqrt(2)
            if sigr > 0:
                xarg = (numpy.range(2*nmk+1)-nmk) / sigr   # exponential arg
                mrkern = numpy.exp(max((-0.5*(xarg**2)),-20.0))
                mrkern = mrkern/mrkern.sum()
            else:
                mrkern = numpy.zeros(2*nmk+1)
                mrkern[nmk] = 1.0    #delta function

    # Construct tangential kernel w/ sigma of sqrt(1-mu**2)*VRT/sqrt(2.)
            if sigt > 0:
                xarg = (numpy.range(2*nmk+1)-nmk) /sigt
                mtkern = exp(max((-0.5*(xarg**2)), -20.0))
                mtkern = mtkern/mtkern.sum()
            else:
                mtkern = numpy.zeros(2*nmk+1)
                mtkern[nmk] = 1.0

    # Sum the radial and tangential components, weighted by surface area
            area_r = 0.5
            area_t = 0.5
            mkern = area_r*mrkern + area_t*mtkern

    # Convolve the total flux profiles, again padding the spectrum on both ends 
    # to protect against Fourier rinnging.
            yfine = scipy.convolve(yfine, mkern, mode='same')
            cfine = scipy.convolve(cfine, mkern, mode='same')

    # Add contribution from current annulus to the running total
            flux += w*yfine
            continuum += w*cfine

        return flux/continuum
