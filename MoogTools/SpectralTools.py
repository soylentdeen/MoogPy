import scipy.signal
import scipy.interpolate
import scipy.optimize
import scipy.integrate
import numpy
import pyfits
import string
import glob

def fitSawtooth(y, window_len=50):
    """
    This routine tries to remove a sawtooth wave from the spectrum.

    As of yet, it DOES NOT work.
    """

    s=numpy.r_[y[window_len-1:0:-1],y,y[-1:-window_len:-1]]
    w = numpy.hanning(window_len)
    newy=numpy.convolve(w/w.sum(),s,mode='valid')
    newy=newy[(window_len/2-1):-(window_len/2)]
    goodPoints = y > newy

    fitfunc = lambda p, x : p[0]+p[1]*scipy.signal.sawtooth(2.*numpy.pi*p[2]*x + p[3])
    errfunc = lambda p, x, y: numpy.abs(fitfunc(p,x) - y)[goodPoints]
    coeffs = [numpy.mean(y), 0.04, 0.0002, 0.0]
    x = numpy.arange(len(newy))
    pfit, success = scipy.optimize.leastsq(errfunc, coeffs, args=(x,newy) )
    print pfit, success
    #raw_input()
    retval = pfit[0]+pfit[1]*scipy.signal.sawtooth(2.*numpy.pi*pfit[2]*x + pfit[3])
    #retval = 0.96-0.04*scipy.signal.sawtooth(2.*numpy.pi*0.0004*x + 0)
    return (retval, newy, goodPoints)
    

def resample(x, y, R, nyquist=False, halt=False):
    """
    This routine convolves a given spectrum to a resolution R
    :INPUTS:
        x: numpy array containing wavelengths
        y: numpy array containing fluxes
        R: Desired resolving power
        nyquist (optional): returns a nyquist-sampled spectrum
               (default: False)

    :RETURNS:
        new_x: new wavelength array
        new_y: new flux array

    :EXAMPLE:
     ::
        highresx, higresy = SpectralTools.read_2col_spectrum('highres.dat')
        lowresx, lowresy = resample(x, y, 2000)
    """
    subsample = 16.0

    xstart = x[0]
    xstop = x[-1]

    newx = [xstart]
    while newx[-1] < xstop:
        stepsize = newx[-1]/(R*subsample)
        newx.append(newx[-1]+stepsize)

    #xdouble = [xstart]
    #while xdouble[-1] < xstop:

    f = scipy.interpolate.interpolate.interp1d(x, y, bounds_error=False)
    newy = f(newx)
    const = numpy.ones(len(newx))

    xk = numpy.array(range(int(4.0*subsample)))
    yk = numpy.exp(-(xk-(2.0*subsample))**2.0/(subsample**2.0/(4.0*numpy.log(2.0))))
    
    result = scipy.signal.convolve(newy, yk, mode ='valid')
    normal = scipy.signal.convolve(const, yk, mode = 'valid')

    bm = numpy.isfinite(result)
    xvals = numpy.array(newx[int(len(xk)/2.0):-int(len(xk)/2.0)])
    yvals = numpy.array(result[bm]/normal[bm])
    if nyquist:
        nyquistx = []
        delta_x = min(xvals)/(2.0*R)
        nyquistx.append(min(xvals)+delta_x)
        while nyquistx[-1] < max(xvals)-delta_x:
            delta_x = nyquistx[-1]/(2.0*R)
            nyquistx.append(nyquistx[-1]+delta_x)
        nyquistx = numpy.array(nyquistx)
        nyquisty = binSpectrum(yvals, xvals, nyquistx)
        retval = (nyquistx, nyquisty)
    else:
        retval = (xvals, yvals)

    if halt:
        print len(xvals), len(yvals), len(result)
        raw_input()
    return retval

def binSpectrum(spectrum, native_wl, new_wl):
    """
        This routine pixelates a synthetic spectrum, in effect simulating the 
        discrete nature of detector pixels.
    """
    
    retval = numpy.zeros(len(new_wl))
    for i in range(len(new_wl)-1):
        bm = scipy.where( (native_wl > new_wl[i]) & (
            native_wl <= new_wl[i+1]))[0]
        if (len(bm) > 1):
            num=scipy.integrate.simps(spectrum[bm], x=native_wl[bm])
            denom = max(native_wl[bm]) - min(native_wl[bm])
            retval[i] = num/denom
        elif (len(bm) == 1):
            retval[i] = 0.0#native_wl[bm]
        else:
            retval[i] = 0.0#retval[-1]
    bm = scipy.where(native_wl > new_wl[-1])[0]
    if len(bm) > 1:
        num = scipy.integrate.simps(spectrum[bm], x=native_wl[bm])
        denom = max(native_wl[bm]) - min(native_wl[bm])
        retval[-1] = num/denom
    else:
        if len(bm) == 1:
            retval[-1] = spectrum[bm]
        else:
            retval[-1] = spectrum[-1]
    return retval

def diff_spectra(x1, y1, x2, y2, pad=False):
    '''
    '''
    x1 = numpy.array(x1)
    y1 = numpy.array(y1)
    x2 = numpy.array(x2)
    y2 = numpy.array(y2)
    overlap_start = numpy.max([numpy.min(x1), numpy.min(x2)])
    overlap_stop = numpy.min([numpy.max(x1), numpy.max(x2)])
    overlap = scipy.where((x1 >= overlap_start) & (x1 <= overlap_stop))

    y = scipy.interpolate.splrep(x2, y2)

    if pad:
        retval = numpy.zeros(len(x1))
        retval[overlap] = y1[overlap] - scipy.interpolate.splev(x1[overlap],y)
        return x1, retval
    else:
        return numpy.array(x1)[overlap], numpy.array(y1)[overlap] - scipy.interpolate.splev(x1[overlap],y)
    
def merge_spectra(x1, y1, x2, y2):
    x1 = numpy.array(x1)
    y1 = numpy.array(y1)
    x2 = numpy.array(x2)
    y2 = numpy.array(y2)
    y1[numpy.isnan(y1)] = 0.0
    y2[numpy.isnan(y2)] = 0.0
    if len(x1) == 0:
        return x2, y2
    if len(x2) == 0:
        return x1, y1
    overlap_start = numpy.max([numpy.min(x1), numpy.min(x2)])
    overlap_stop = numpy.min([numpy.max(x1), numpy.max(x2)])
    overlap = scipy.where((x1 >= overlap_start) & (x1 <= overlap_stop))
    if (len(overlap[0]) > 1):
        unique1 = scipy.where((x1 < overlap_start) | (x1 > overlap_stop))
        unique2 = scipy.where((x2 < overlap_start) | (x2 > overlap_stop))

        y = scipy.interpolate.splrep(x2, y2)

        new_x = numpy.append(x1, x2[unique2])
        yinterp = scipy.interpolate.splev(x1[overlap], y)
        merged = y1[overlap] + yinterp
        new_y = numpy.append(numpy.append(y1[unique1], merged), y2[unique2])
    else:
        new_x = numpy.append(x1, x2)
        new_y = numpy.append(y1, y2)

    return new_x, new_y


def interpolate_spectrum(x1, x2, y1, pad=None):
    overlap_start = numpy.max([numpy.min(x1), numpy.min(x2)])
    overlap_stop = numpy.min([numpy.max(x1), numpy.max(x2)])
    overlap = ( x2 >= overlap_start) & (x2 <= overlap_stop)

    y = scipy.interpolate.splrep(x1, y1)
    if pad!=None:
        retval = numpy.ones(len(x2))* pad
        retval[overlap] = scipy.interpolate.splev(x2[overlap],y)
        return retval
    else:
        return scipy.interpolate.splev(x2[overlap],y)

def write_2col_spectrum(filename, wl, fl):
    '''
    Prints a spectrum to a two-column data file
    '''

    data = open(filename, 'w')

    for line in zip(wl, fl):
        data.write(str(line[0])+' '+str(line[1])+'\n')

    data.close()

def write_3col_spectrum(filename, wl, fl, er):
    '''
    Prints a spectrum to a three-column data file
    '''

    data = open(filename, 'w')

    for line in zip(wl, fl, er):
        data.write(str(line[0])+' '+str(line[1])+' '+str(line[2])+'\n')

    data.close()

def write_MOOG_obs_spectrum(filename, wl, fl):
    '''
    Prints a spectrum to a two-column data file in a format useful for the MOOG spectral synthesis program.
    '''

    data = open(filename, 'w')

    for d in zip(wl, fl):
        data.write('%10f %9f\n' % (d[0], d[1]) )

    data.close()

def read_2col_spectrum(filename):
    '''
    Reads in a spectrum from a text file in a 2 column format.

    column 1: wavelength
    column 2: flux
    '''
    data = open(filename).read().split('\n')
    
    x = []
    y = []
    
    for line in data:
        l = line.split()
        if len(l) == 2:
            x.append(float(l[0]))
            y.append(float(l[1]))
            
    x = numpy.array(x)
    y = numpy.array(y)
    
    return x, y

def read_3col_spectrum(filename):
    '''
    Reads in a spectrum from a text file in a 2 column format.

    column 1: wavelength
    column 2: flux
    '''
    data = open(filename).read().split('\n')
    
    x = []
    y = []
    z = []
    
    for line in data:
        l = line.split()
        if len(l) == 3:
            x.append(float(l[0]))
            y.append(float(l[1]))
            z.append(float(l[2]))
            
    x = numpy.array(x)
    y = numpy.array(y)
    z = numpy.array(z)
    
    return x, y, z

class spectrumUnit( object ):
    def __init__(self, wl=None, I=None, Q=None, U=None, V=None):
        self.wl = wl
        self.flux_I = I
        self.flux_Q = Q
        self.flux_U = U
        self.flux_V = V
        
        
    def bin(self, newWl):
        """
            This routine simulates the discrete nature of detector pixels.
        """
    
        if self.flux_I != None:
            newSpec_I = numpy.zeros(len(newWl))
        if self.flux_Q != None:
            newSpec_Q = numpy.zeros(len(newWl))
        if self.flux_U != None:
            newSpec_U = numpy.zeros(len(newWl))
        if self.flux_V != None:
            newSpec_V = numpy.zeros(len(newWl))
        for i in range(len(newWl)-1):
            inBin = scipy.where( (self.wl > newWl[i]) & (
                self.wl <= newWl[i+1]))[0]
            if (len(inBin) > 1):
                denom = self.wl[inBin][-1] - self.wl[inBin][0]
                if self.flux_I != None:
                    num=scipy.integrate.simps(self.flux_I[inBin], 
                            x=self.wl[inBin])
                    newSpec_I[i] = num/denom
                if self.flux_Q != None:
                    num=scipy.integrate.simps(self.flux_Q[inBin], 
                            x=self.wl[inBin])
                    newSpec_Q[i] = num/denom
                if self.flux_U != None:
                    num=scipy.integrate.simps(self.flux_U[inBin], 
                            x=self.wl[inBin])
                    newSpec_U[i] = num/denom
                if self.flux_V != None:
                    num=scipy.integrate.simps(self.flux_V[inBin], 
                            x=self.wl[inBin])
                    newSpec_V[i] = num/denom
            elif (len(inBin) == 1):
                if self.flux_I != None:
                    newSpec_I[i] = 0.0
                if self.flux_Q != None:
                    newSpec_Q[i] = 0.0
                if self.flux_U != None:
                    newSpec_U[i] = 0.0
                if self.flux_V != None:
                    newSpec_V[i] = 0.0
            else:
                if self.flux_I != None:
                    newSpec_I[i] = 0.0
                if self.flux_Q != None:
                    newSpec_Q[i] = 0.0
                if self.flux_U != None:
                    newSpec_U[i] = 0.0
                if self.flux_V != None:
                    newSpec_V[i] = 0.0

        inBin = scipy.where(self.wl > newWl[-1])[0]
        if len(inBin) > 1:
            denom = self.wl[inBin][-1] - self.wl[inBin][0]
            if self.flux_I != None:
                num = scipy.integrate.simps(self.flux_I[inBin], 
                        x=self.wl[inBin])
                newSpec_I[-1] = num/denom
            if self.flux_Q != None:
                num = scipy.integrate.simps(self.flux_Q[inBin], 
                        x=self.wl[inBin])
                newSpec_Q[-1] = num/denom
            if self.flux_U != None:
                num = scipy.integrate.simps(self.flux_U[inBin], 
                        x=self.wl[inBin])
                newSpec_U[-1] = num/denom
            if self.flux_V != None:
                num = scipy.integrate.simps(self.flux_V[inBin], 
                        x=self.wl[inBin])
                newSpec_V[-1] = num/denom
        else:
            if len(inBin) == 1:
                if self.flux_I != None:
                    newSpec_I[-1] = self.flux_I[inBin]
                if self.flux_Q != None:
                    newSpec_Q[-1] = self.flux_Q[inBin]
                if self.flux_U != None:
                    newSpec_U[-1] = self.flux_U[inBin]
                if self.flux_V != None:
                    newSpec_V[-1] = self.flux_V[inBin]
            else:
                if self.flux_I != None:
                    newSpec_I[-1] = self.flux_I[-1]
                if self.flux_Q != None:
                    newSpec_Q[-1] = self.flux_Q[-1]
                if self.flux_U != None:
                    newSpec_U[-1] = self.flux_U[-1]
                if self.flux_V != None:
                    newSpec_V[-1] = self.flux_V[-1]

        self.wl = newWl
        if self.flux_I != None:
            self.flux_I = newSpec_I
        if self.flux_Q != None:
            self.flux_Q = newSpec_Q
        if self.flux_U != None:
            self.flux_U = newSpec_U
        if self.flux_V != None:
            self.flux_V = newSpec_V

    def __sub__(self, other):
        return None

class Spectrum( object ):
    def __init__(self, name):
        self.name = None
        self.native = spectrumUnit()
        self.processed = spectrumUnit()
        self.wlStart = None
        self.wlStop = None

    def resample(self, R, nyquist=False):
        """
        This routine convolves a given spectrum to a resolution R
        :INPUTS:
            R: Desired resolving power
            nyquist (optional): returns a nyquist-sampled spectrum
                   (default: False)
    
        :RETURNS:
            new_x: new wavelength array
            new_y: new flux array

        :EXAMPLE:
         ::
            highres = SpectralTools.Spectrum('highres.dat')
            highres.resample(2000)
        """
        subsample = 16.0

        newWl = [self.wlStart]
        while newWl[-1] < self.wlStop:
            stepsize = newWl[-1]/(R*subsample)
            newWl.append(newWl[-1]+stepsize)

        if self.native.flux_I != None:
            I = scipy.interpolate.interpolate.interp1d(self.native.wl,
                    self.native.flux_I, bounds_error=False)
            newI = I(newWl)
        if self.native.flux_V != None:
            V = scipy.interpolate.interpolate.interp1d(self.native.wl,
                    self.native.flux_V bounds_error=False)
            newV = V(newWl)
        if self.native.flux_Q != None:
            Q = scipy.interpolate.interpolate.interp1d(self.native.wl,
                    self.native.flux_Q bounds_error=False)
            newQ = Q(newWl)
        if self.native.flux_U != None:
            U = scipy.interpolate.interpolate.interp1d(self.native.wl,
                    self.native.flux_U bounds_error=False)
            newU = U(newWl)
        const = numpy.ones(len(newWl))

        xk = numpy.array(range(int(4.0*subsample)))
        yk = numpy.exp(-(xk-(2.0*subsample))**2.0/(subsample**2.0/(4.0*numpy.log(2.0))))
    
        newWl = numpy.array(newWl[int(len(xk)/2.0):-int(len(xk)/2.0)])

        normal = scipy.signal.convolve(const, yk, mode = 'valid')
        result_I = scipy.signal.convolve(newI, yk, mode ='valid')
        goodPoints = numpy.isfinite(result_I)
        flux_I = numpy.array(result_I[goodPoints]/normal[goodPoints])
        if self.native_flux_V != None:
            result_V = scipy.signal.convolve(newV, yk, mode ='valid')
            flux_V = numpy.array(result_V[goodPoints]/normal[goodPoints])
        else:
            flux_V = None
        if self.native_flux_Q != None:
            result_Q = scipy.signal.convolve(newQ, yk, mode ='valid')
            flux_Q = numpy.array(result_Q[goodPoints]/normal[goodPoints])
        else:
            flux_Q = None
        if self.native_flux_U != None:
            result_U = scipy.signal.convolve(newU, yk, mode ='valid')
            flux_U = numpy.array(result_U[goodPoints]/normal[goodPoints])
        else:
            flux_U = None

        self.processed = SpectrumUnit(newWl, flux_I, flux_Q, flux_U, flux_V)

        if nyquist:
            nyquistWl = []
            deltaWl = min(newWl)/(2.0*R)
            nyquistWl.append(min(newWl) + deltaWl)
            while nyquistWl[-1] < max(newWl) - deltaWl:
                deltaWl = nyquistWl[-1]/(2.0*R)
                nyquistWl.append(nyquistWl[-1] + deltaWl)
            newWl = numpy.array(nyquistWl)
            self.processed.binSpectrum(newWl=newWl)
    


class MoogStokesSpectrum( Spectrum ):
    def __init__(self, name='', memory=False, **kwargs):
        self.memory = memory
        if self.memory:
            self.parent = kwargs["PARENT"]
            self.deltav = self.parent.deltav
            self.vsini = self.parent.vsini
        else:
            self.deltav = kwargs["DELTAV"]
            self.vsini = kwargs["VSINI"]

            self.name = name
            self.info = pyfits.info(self.name, output='')
            self.nheaders = len(self.info)-1
            self.primaryHeader = pyfits.getheader(f, ext=0)
            self.headers = []
            self.wavelengthRanges = []
            for i in range(self.nheaders):
                self.headers.append(pyfits.getheader(self.name, ext=i+1))
                self.wavelengthRanges.append([self.headers[-1].get('WLSTART'),
                    self.headers[-1].get('WLSTOP')])

    def integrate(self, **kwargs):
        self.loadAngles()
        self.loadSpectra(kwargs)
        self.interpolateSpectra()
        self.diskInt()

    def loadAngles(self):
        if self.memory:
            self.phi = self.parent.phi_angle[:self.parent.ncells]
            self.mu = self.parent.mus[:self.parent.ncells]
        else:
            self.phi = []
            self.mu = []
            for hdr in self.headers:
                self.phi.append(hdr.get('PHI_ANGLE'))
                self.mu.append(hdr.get('MU'))

    def loadSpectra(self, **kwargs):
        if self.memory:
            self.I = numpy.array(self.parent.flux_I)/numpy.array(self.parent.continuum)
            self.V = numpy.array(self.parent.flux_V)/numpy.array(self.parent.continuum)
            self.continuum = numpy.array(self.parent.continuum)
            self.wl = numpy.array(self.parent.wave)
        else:
            if "wlRange" in kwargs.keys():
                self.wlStart = kwargs["wlRange"][0]
                self.wlStop = kwargs["wlRange"][1]



def winnow_MoogStokes_Spectra(directory, wlStart, wlStop, trackedParams=None):
    files = glob.glob(directory+'*.fits')
    header_keys = trackedParams.keys()
    waves = []
    fluxes = []
    for f in files:
        info = pyfits.info(f, output='')
        nheaders = len(info)-1
        primary_header = pyfits.getheader(f, ext=0)
        for i in range(nheaders):
            hdr = pyfits.getheader(f, ext=i+1)
            if (hdr.get('WLSTART') < wlStop) and (hdr.get('WLSTOP') > wlStart):
                data = pyfits.getdata(f, ext=i+1)
                waves.append(data.field('Wavelength'))
                fluxes.append(data.field('Stokes_I'))
                for hk in header_keys:
                    trackedParams[hk].append(primary_header.get(hk))
    return trackedParams, waves, fluxes


def read_IGRINS_spectrum(starFile, A0VFile):
    star_Hdu = pyfits.open(starFile, ignore_missing_end=True)
    A0V_Hdu = pyfits.open(A0VFile, ignore_missing_end=True)
    
    star_hdr = star_Hdu[0].header
    star = star_Hdu[0].data
    
    A0V_hdr = A0V_Hdu[0].header
    A0V = A0V_Hdu[0].header
    
    
    

def read_fits_spectrum(filename):
    hdulist = pyfits.open(filename, ignore_missing_end=True)
    hdr = hdulist[0].header
    dat = hdulist[0].data

    wl = dat[0]
    fl = dat[1]
    dfl = dat[2]
    
    bm = scipy.where( numpy.isfinite(fl) == True)

    return (wl[bm], fl[bm], dfl[bm])

def read_IRAF_fits_spectrum(filename):
    """ Reads in an echelle spectrum which has been reduced by IRAF """
    hdulist = pyfits.open(filename, ignore_missing_end = True)
    hdr = hdulist[0].header
    dat = hdulist[0].data

    "Finds the number of orders"
    #n_orders = int(hdr["NAXIS2"])
    waveTable = hdr["WAT2*"]

    "Strings together the wavelength conversion strings"
    linear_header = ''
    for tableEntry in waveTable.items():
        linear_header += string.ljust(tableEntry[1], 68)

    """
    Extracts the coefficients necessary for the wavelength solution
    for each order
    """
    wlsol = []
    for sol in linear_header.split('spec')[2:]:
        sol = sol.split()
        wlsol.append([float(sol[5]), float(sol[6])])

    """
    Creates the wavelength solution for each pixel, and places it in a
    parallel array.
    """
    orders = []
    for i in range(len(dat)):
        wl = numpy.arange(len(dat[i]))*wlsol[i][1]+wlsol[i][0]
        orders.append([wl, dat[i]])

    orders = numpy.array(orders)
    return hdr, orders

def fit_gaussians(x, y, linecenters, R, **kwargs):

    params = []
    strength = -0.05
    for i in range(len(linecenters)):
        fwhm = 0.05*linecenters[i]/R
        if "strengthGuesses" in kwargs:
            params.append(kwargs["strengthGuesses"][i])
        else:
            params.append(strength)      #Strength
        params.append(linecenters[i])
        params.append(fwhm)          #FWHM

    def fitfunc(pars, xpts):
        retval = numpy.ones(len(xpts))
        for i in range(len(linecenters)):
            k = i*3
            for j in range(len(xpts)):
                retval[j] += pars[k]*numpy.exp(-(xpts[j]-pars[k+1])**2.0/(pars[k+2]))
        return retval

    def errfunc(pars, xpts, ypts):
        return numpy.abs(fitfunc(pars, xpts) - ypts)
        
    pfit = scipy.optimize.leastsq(errfunc, params, args=(x, y))

    coeffs = pfit[0]

    fit = numpy.ones(len(x))
    for i in range(len(linecenters)):
        k = i*3
        for j in range(len(x)):
            fit[j] += coeffs[k]*numpy.exp(-(x[j]-coeffs[k+1])**2.0/(coeffs[k+2]))

    return coeffs, fit

def calc_EW(x, y, xstart, xstop):
   bm = scipy.where( (x > xstart) & (x < xstop) )[0]
   cont = numpy.ones(len(bm))
   num = scipy.integrate.simps(y[bm], x[bm])
   denom = scipy.integrate.simps(cont, x[bm])
   return (denom-num)


def blackBody(**kwargs):
    """ Returns a blackbody function over the given wavelength  """

    """
    inputs:
        wl : wavelength array
               Assumed to be in units of cm, unless specified by wlUnits kwarg

        nu : frequency array
               
        T : Blackbody Temperature (K)

        wlUnits: Units of wavelengths (Optional)
               'Angstroms'
               'nanometers'
               'microns'
               'cm'
               'meters'
        outUnits: cgs Units to output.
               'Fnu'
               'Flambda'
               'Energy'

    outputs:
        y : Blackbody function.  Unit is assumed to be Flambda or Fnu (depnding on whether
                wl or nu was given) unless overridden by outUnits kwarg
    """
    h = 6.626e-27
    c = 2.998e10
    k = 1.38e-16
    T = kwargs["T"]

    if "wl" in kwargs:
        wl = kwargs["wl"]
        c1 = 2.0*h*c**2.0
        c2 = h*c/(k*T)
        Flambda = c1/(wl**5.0*(numpy.exp(c2/wl)-1.0))
        if "outUnits" in kwargs:
            if kwargs["outUnits"] == "Energy":
                return Flambda*wl
        else:
            return Flambda
    elif "nu" in kwargs:
        nu = kwargs["nu"]
        c1 = 2.0*h/(c**2.0)
        c2 = h/(k*T)
        Fnu = c1*nu**2.0/(numpy.exp(c2*nu) - 1.0)
        if "outUnits" in kwargs:
            if kwargs["outUnits"] == "Energy":
                return Fnu*nu
        else:
            return Fnu


class photometrySynthesizer( object ):
    def __init__(self, **kwargs):
        if "filterDir" in kwargs:
            self.fdir = filter_dir
        else:
            self.fdir = '/home/deen/Data/StarFormation/Photometry/FILTER_PROFILES/'

        filterNames = ['Uj', 'Bj', 'Vj', 'Rc', 'Ic', '2massj', '2massh', '2massk']
        fileNames = ['U_Landolt.dat', 'B_Bessell.dat', 'V_Bessell.dat', 'cousins_Rband.dat', 'cousins_Iband.dat', 'J_2MASS.dat', 'H_2MASS.dat', 'K_2MASS.dat']
        fnu_zero = [1829, 4144, 3544, 2950, 2280.0, 1594.0, 1024.0, 666.7 ]
        flam_zero = [4.0274905e-09, 6.3170333e-09, 3.6186341e-09, 2.1651655e-9, 1.1326593e-09, 3.129e-10, 1.133e-10, 4.283e-11] #erg/s/cm^2/Angstrom
        lambda_eff = [3600, 4362, 5446, 6413, 7978, 12285, 16385, 21521]
        mVega = [0.02, 0.02, 0.03, 0.039, 0.035, -0.001, +0.019, -0.017]

        self.photBands = []
        for band in zip(filterNames, fileNames, fnu_zero, flam_zero, lambda_eff, mVega):
            photBand = dict()
            photBand['Name'] = band[0]
            photBand['file'] = band[1]
            photBand['fnu_zero'] = band[2]
            photBand['flam_zero'] = band[3]
            photBand['lambda_eff'] = band[4]
            photBand['mVega'] = band[5]
            
            fx = []
            fy = []
            dat = open(self.fdir+band[1], 'r').read().split('\n')
            for line in dat:
                if len(line) > 0:
                    l = line.split()
                    fx.append(float(l[0])*1e4)
                    fy.append(float(l[1]))
            fy = numpy.array(fy)
            fx = numpy.array(fx)
            photBand['min_x'] = min(fx)
            photBand['max_x'] = max(fx)
            photBand['photSpline'] = scipy.interpolate.splrep(fx, fy)

            self.photBands.append(photBand)
    
    def photFlux(self, x, y, filtName):
        for band in self.photBands:
            if band['Name'] == filtName:
                bm = scipy.where( (x > band['min_x'] ) & (x < band['max_x']) )[0]
                fnew = scipy.interpolate.splev(x[bm], band['photSpline'])
                valid_bm = scipy.where( (fnew > 0.0) & (y[bm] > 0.0) )[0]
                numerator = scipy.integrate.simps(y[bm][valid_bm]*fnew[valid_bm], x[bm][valid_bm])
                denom = scipy.integrate.simps(fnew[valid_bm], x[bm][valid_bm])
                return numerator/denom
