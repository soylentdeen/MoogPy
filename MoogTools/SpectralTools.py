import scipy.signal
import scipy.interpolate
import scipy.optimize
import scipy.integrate
import numpy
import astropy.io.fits as pyfits
import string
import glob
import os
import string
import random

class SpectrumError( Exception ):
    def __init__(self, value, errmsg):
        '''
        SpectrumError
        
        Raised on execptions within Spectrum objects
        
        Error definitions are as follows:
        0 | Failure loading Raw Data!!
        1 | Failure loading Processed Data!!
        2 | Failure calculating Equivalent Widths!
        3 | Failure Convolving Spectrum!
        4 | Failure Calculating Difference Spectrum!
        '''
        self.value = value
        self.message = {}
        self.message[0] = "Failure loading Raw Data!! %s" % errmsg
        self.message[1] = "Failure loading Processed Data!! %s" % errmsg
        self.message[2] = "Failure calculating Equivalent Widths! %s" % errmsg
        self.message[3] = "Failure Convolving Spectrum! %s" % errmsg
        self.message[4] = "Failure Calculating Difference Spectrum! %s" % errmsg

    def __str__(self):
        return repr(self.message[self.value])

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
    #print pfit, success
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
        print("%d %d %d" % (len(xvals), len(yvals), len(result)))
        input()
    return retval

def binSpectrum(spectrum, new_wl):

    overlap_start = numpy.max([numpy.min(spectrum.wl), numpy.min(new_wl)])
    overlap_stop = numpy.min([numpy.max(spectrum.wl), numpy.max(new_wl)])
    overlap = scipy.where((new_wl >= overlap_start) & (new_wl <= overlap_stop))
    
    new_wl = new_wl[overlap]
    
    if spectrum.flux_I != None:
        new_I = numpy.zeros(len(new_wl))
    else:
        new_I = None
    if spectrum.flux_Q != None:
        new_Q = numpy.zeros(len(new_wl))
    else:
        new_Q = None
    if spectrum.flux_U != None:
        new_U = numpy.zeros(len(new_wl))
    else:
        new_U = None
    if spectrum.flux_V != None:
        new_V = numpy.zeros(len(new_wl))
    else:
        new_V = None
    if spectrum.continuum != None:
        new_continuum = numpy.zeros(len(new_wl))
    else:
        new_continuum = None

    for i in range(len(new_wl)-1):
        bm = scipy.where( (spectrum.wl > new_wl[i]) & (
                spectrum.wl <= new_wl[i+1]))[0]
        if (len(bm) > 1):
            denom = max(spectrum.wl[bm]) - min(spectrum.wl[bm])
            if spectrum.flux_I != None:
                num = scipy.integrate.simps(spectrum.flux_I[bm], x=spectrum.wl[bm])
                new_I[i] = num/denom
            if spectrum.flux_Q != None:
                num = scipy.integrate.simps(spectrum.flux_Q[bm], x=spectrum.wl[bm])
                new_Q[i] = num/denom
            if spectrum.flux_U != None:
                num = scipy.integrate.simps(spectrum.flux_U[bm], x=spectrum.wl[bm])
                new_U[i] = num/denom
            if spectrum.flux_V != None:
                num = scipy.integrate.simps(spectrum.flux_V[bm], x=spectrum.wl[bm])
                new_V[i] = num/denom
            if spectrum.continuum != None:
                num = scipy.integrate.simps(spectrum.continuum[bm], x=spectrum.wl[bm])
                new_continuum[i] = num/denom
        elif (len(bm) == 1):
            if spectrum.flux_I != None:
                new_I[i] = spectrum.flux_I[bm]
            if spectrum.flux_Q != None:
                new_Q[i] = spectrum.flux_Q[bm]
            if spectrum.flux_U != None:
                new_U[i] = spectrum.flux_U[bm]
            if spectrum.flux_V != None:
                new_V[i] = spectrum.flux_V[bm]
            if spectrum.continuum != None:
                new_continuum[i] = spectrum.continuum[bm]
        else:
            if spectrum.flux_I != None:
                new_I[i] = 0.0 
            if spectrum.flux_Q != None:
                new_Q[i] = 0.0 
            if spectrum.flux_U != None:
                new_U[i] = 0.0 
            if spectrum.flux_V != None:
                new_V[i] = 0.0 
            if spectrum.continuum != None:
                new_continuum[i] = 0.0 
    bm = scipy.where(spectrum.wl > new_wl[-1])[0]
    if len(bm) > 1:
        denom = max(spectrum.wl[bm]) - min(spectrum.wl[bm])
        if spectrum.flux_I != None:
            num = scipy.integrate.simps(spectrum.flux_I[bm], x=spectrum.wl[bm])
            new_I[-1] = num/denom
        if spectrum.flux_Q != None:
            num = scipy.integrate.simps(spectrum.flux_Q[bm], x=spectrum.wl[bm])
            new_Q[-1] = num/denom
        if spectrum.flux_U != None:
            num = scipy.integrate.simps(spectrum.flux_U[bm], x=spectrum.wl[bm])
            new_U[-1] = num/denom
        if spectrum.flux_V != None:
            num = scipy.integrate.simps(spectrum.flux_V[bm], x=spectrum.wl[bm])
            new_V[-1] = num/denom
        if spectrum.continuum != None:
            num = scipy.integrate.simps(spectrum.continuum[bm], x=spectrum.wl[bm])
            new_continuum[-1] = num/denom
    else:
        if len(bm) == 1:
            if spectrum.flux_I != None:
                new_I[-1] = spectrum.flux_I[bm]
            if spectrum.flux_Q != None:
                new_Q[-1] = spectrum.flux_Q[bm]
            if spectrum.flux_U != None:
                new_U[-1] = spectrum.flux_U[bm]
            if spectrum.flux_V != None:
                new_V[-1] = spectrum.flux_V[bm]
            if spectrum.continuum != None:
                new_continuum[-1] = spectrum.continuum[bm]
        else:
            if spectrum.flux_I != None:
                new_I[-1] = spectrum.flux_I[-1]
            if spectrum.flux_Q != None:
                new_Q[-1] = spectrum.flux_Q[-1]
            if spectrum.flux_U != None:
                new_U[-1] = spectrum.flux_U[-1]
            if spectrum.flux_V != None:
                new_V[-1] = spectrum.flux_V[-1]
            if spectrum.continuum != None:
                new_continuum[-1] = spectrum.continuum[-1]
                
    retval = Spectrum(wl=new_wl, I=new_I, Q=new_Q,
                U=new_U, V=new_V, continuum=new_continuum, 
                header=spectrum.header, spectrum_type="REBINNED",
                filename=spectrum.filename, ext=spectrum.ext, 
                preserve=spectrum.preserve)
                
    return retval
        
def binSpectrum_old(spectrum, native_wl, new_wl):
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
    
def mergeSpectra(first=None, second=None):
    if second == None:
        return first
    if first == None:
        return second
    
    x1 = first.wl
    x2 = second.wl

    overlap_start = numpy.max([numpy.min(x1), numpy.min(x2)])
    overlap_stop = numpy.min([numpy.max(x1), numpy.max(x2)])
    overlap = scipy.where((x1 >= overlap_start) & (x1 <= overlap_stop))
    if (len(overlap[0]) > 1):
        unique1 = scipy.where((x1 < overlap_start) | (x1 > overlap_stop))
        unique2 = scipy.where((x2 < overlap_start) | (x2 > overlap_stop))

        new_x = numpy.append(x1, x2[unique2])

        if (first.flux_I != None) & (second.flux_I != None):
            I1 = first.flux_I
            I2 = second.flux_I
            I1[numpy.isnan(I1)] = 0.0
            I2[numpy.isnan(I2)] = 0.0
            I = scipy.interpolate.splrep(x2, I2)
            Iinterp = scipy.interpolate.splev(x1[overlap], I)
            mergedI = I1[overlap] + Iinterp
            new_I = numpy.append(numpy.append(I1[unique1], merged), I2[unique2])
        else:
            new_I = None
        if (first.flux_Q != None) & (second.flux_Q != None):
            Q1 = first.flux_Q
            Q2 = second.flux_Q
            Q1[numpy.isnan(Q1)] = 0.0
            Q2[numpy.isnan(Q2)] = 0.0
            Q = scipy.interpolate.splrep(x2, Q2)
            Qinterp = scipy.interpolate.splev(x1[overlap], Q)
            mergedQ = Q1[overlap] + Qinterp
            new_Q = numpy.append(numpy.append(Q1[unique1], merged), Q2[unique2])
        else:
            new_Q = None
        if (first.flux_U != None) & (second.flux_U != None):
            U1 = first.flux_U
            U2 = second.flux_U
            U1[numpy.isnan(U1)] = 0.0
            U2[numpy.isnan(U2)] = 0.0
            U = scipy.interpolate.splrep(x2, U2)
            Uinterp = scipy.interpolate.splev(x1[overlap], U)
            mergedU = U1[overlap] + Uinterp
            new_U = numpy.append(numpy.append(U1[unique1], merged), U2[unique2])
        else:
            new_U = None
        if (first.flux_V != None) & (second.flux_V != None):
            V1 = first.flux_V
            V2 = second.flux_V
            V1[numpy.isnan(V1)] = 0.0
            V2[numpy.isnan(V2)] = 0.0
            V = scipy.interpolate.splrep(x2, V2)
            Vinterp = scipy.interpolate.splev(x1[overlap], V)
            mergedV = V1[overlap] + Vinterp
            new_V = numpy.append(numpy.append(V1[unique1], merged), V2[unique2])
        else:
            new_V = None
        if (first.continuum != None) & (second.continuum != None):
            C1 = first.continuum
            C2 = second.continuum
            C1[numpy.isnan(C1)] = 0.0
            C2[numpy.isnan(C2)] = 0.0
            C = scipy.interpolate.splrep(x2, C2)
            Cinterp = scipy.interpolate.splev(x1[overlap], C)
            mergedC = C1[overlap] + Cinterp
            new_C = numpy.append(numpy.append(C1[unique1], merged), C2[unique2])
        else:
            new_C = None
            if new_I != None:
                new_I /= 2.0
            if new_Q != None:
                new_Q /= 2.0
            if new_U != None:
                new_U /= 2.0
            if new_V != None:
                new_V /= 2.0
    else:
        new_x = numpy.append(x1, x2)
        if (first.flux_I != None) & (second.flux_I != None):
            I1 = first.flux_I
            I2 = second.flux_I
            I1[numpy.isnan(I1)] = 0.0
            I2[numpy.isnan(I2)] = 0.0
            new_I = numpy.append(I1, I2)
        else:
            new_I = None
        if (first.flux_Q != None) & (second.flux_Q != None):
            Q1 = first.flux_Q
            Q2 = second.flux_Q
            Q1[numpy.isnan(Q1)] = 0.0
            Q2[numpy.isnan(Q2)] = 0.0
            new_Q = numpy.append(Q1, Q2)
        else:
            new_Q = None
        if (first.flux_U != None) & (second.flux_U != None):
            U1 = first.flux_U
            U2 = second.flux_U
            U1[numpy.isnan(U1)] = 0.0
            U2[numpy.isnan(U2)] = 0.0
            new_U = numpy.append(U1, U2)
        else:
            new_U = None
        if (first.flux_V != None) & (second.flux_V != None):
            V1 = first.flux_V
            V2 = second.flux_V
            V1[numpy.isnan(V1)] = 0.0
            V2[numpy.isnan(V2)] = 0.0
            new_V = numpy.append(V1, V2)
        else:
            new_V = None
        if (first.continuum != None) & (second.continuum != None):
            C1 = first.continuum
            C2 = second.continuum
            C1[numpy.isnan(C1)] = 0.0
            C2[numpy.isnan(C2)] = 0.0
            new_C = numpy.append(C1, C2)
        else:
            new_C = None
        
    retval = Spectrum(wl=new_x, I = new_I, Q = new_Q, U = new_U, V = new_V, 
            continuum = new_C, spectrum_type='MERGED')

    return retval

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

class Spectrum( object ):
    def __init__(self, wl=None, I=None, dI=None, Q=None, U=None, V=None,
            continuum=None, header=pyfits.Header(), spectrum_type=None,
            filename=None, ext=None, preserve=False, label=None):
        """
            Spectrum.__init__(wl=None, I=None, Q=None, U=None, V=None,
            continuum=None, header=pyfits.Header(), spectrum_type=None,
            filename=None, ext=None, preserve=False)
            
            Creates a Spectrum object from arrays
            
            wl = numpy array of wavelength values (should be microns?)
            I = numpy array of Stokes I flux values
            Q = numpy array of Stokes Q flux values
            U = numpy array of Stokes U flux values
            V = numpy array of Stokes V flux values
            continuum = numpy array of continuum values.  Can/should be omitted
                if the flux values are normalized.
            header = pyfits Header object containing the FITS header to be 
                saved with the Spectrum when it is saved
            spectrum_type = 'RAW', 'INTERPOLATED', 'DISK INTEGRATED',
                'CONVOLVED', 'BLENDED', 'MERGED', 'ROTATED', 'DIFFERENCE',
                'SCALED', 'MOOG DISK INTEGRATED', 'MOOG EMERGENT'
            filename = Filename for saving (or reading)
            ext = For FITS files with multiple parts, ext is the extention number
            preserve = When True, prepares the fluxes and wavelengths for saving into
                a FITS table.
        """
        self.wl = wl
        self.flux_I = I
        self.dflux_I = dI
        self.flux_Q = Q
        self.flux_U = U
        self.flux_V = V
        self.continuum = continuum
        self.header = header
        self.filename = filename
        self.ext = ext
        self.label = label
        if spectrum_type != None:
            self.addHistory(spectrum_type=spectrum_type)
        if preserve == True:
            self.preserve(I=I!=None, dI=dI!=None, Q=Q!=None, U=U!=None, V=V!=None, continuum=continuum!=None)

    @classmethod
    def from_file(self, header=None, data=None, filename=None, ext=None, label=None):
        """
            Spectrum.from_file(header=None, data=None, filename=None, ext=None)
            
            Creates a Spectrum object by reading previously saved data from a file
            
            header = pyfits Header object to be used in place of any existing FITS
                header.
            data = object containing pyfits datafields (i.e. pyfits.getdata())
            filename = full path of file containing data
            ext = if multiple spectra are contained in the same file, ext 
                specifies which extension to grab.
        """
        if not(data==None):
            self.extractData(data)
        else:
            wl = None
            I = None
            dI = None
            Q = None
            U = None
            V = None
            continuum=None
            if (filename==None) and (ext==None):
                errmsg = "Filename or extension not provided!"
                raise SpectrumError(0, errmsg)

        return self(wl=wl, I=I, dI=dI, Q=Q, U=U, V=V, continuum=continuum, header=header,
                filename=filename, ext=ext, label=label)

    def addHistory(self, spectrum_type=""):
        """
            Spectrum.addHistory(spectrum_type="")
            
            addHistory allows the user to keep track of the provenance of the data
                contained inside.  If an operation changes the type of data (say by 
                convolving to a different resolution), this routine stores the
                "parent's" SPECTRUM ID in the history.  Then, the spectrum is given
                a new SPECTRUM_ID by generating a random string of ascii characters.
        """
        if 'SPECTRUM_ID' in self.header.iterkeys():
            self.header.add_history(self.header.get('SPECTRUM_TYPE')+
                    ' - '+self.header.get('SPECTRUM_ID'))
        self.header.set('SPECTRUM_ID', ''.join(random.choice(string.ascii_letters)
            for _ in range(10)))
        self.header.set('SPECTRUM_TYPE', spectrum_type)
    
    def addLabel(self, label=None):
        self.label = label

    def extractData(self, data):
        """
        Spectrum.extractData(data)
        
        This routine extracts spectrum data from a FITS binary table.
        
        data = pyfits binary table
        
        Note: the data object must contain a 'Wavelength' field.
        """
        try:
            self.wl = data.field('Wavelength')
        except:
            raise SpectrumError(0, "Spectrum data must contain WL data")
        try:
            self.flux_I = data.field('Stokes_I')
        except:
            self.flux_I = None
        try:
            self.dflux_I = data.field('dStokes_I')
        except:
            self.dflux_I = None
        try:
            self.flux_Q = data.field('Stokes_Q')
        except:
            self.flux_Q = None
        try:
            self.flux_U = data.field('Stokes_U')
        except:
            self.flux_U = None
        try:
            self.flux_V = data.field('Stokes_V')
        except:
            self.flux_V = None
        try:
            self.continuum = data.field('Continuum')
        except:
            self.continuum = None

    def loadData(self):
        """
        Spectrum.loadData()
        
        Loads the data related to a spectrum stored in the self.filename.
        
        Once the data is loaded from the file, the file is closed, and 
        removed from memory.  Then, the data is extracted and loaded into the
        Spectrum object
        """
        try:
            datafile = open(self.filename, 'rb')
            data = pyfits.getdata(datafile, ext=self.ext, memmap=False)
            datafile.close()
            del(datafile)
        except:
            raise SpectrumError(0, "Error reading extension %d from %s" %
                    (self.ext, self.filename))
        self.extractData(data)
        del(data)

    def preserve(self, prepareColumns=True, I=True, dI=False, Q=False, U=False, V=True, continuum=True):
        """
        Spectrum.preserve(prepareColumns=True, I=True, Q=False, U=False, V=True, continuum=True)
        
        prepares the spectrum for saving in a FITS binary table by
        creating the columns object (a pyfits.ColDefs object) from
        the Spectrum data.
        
        prepareColumns [Boolean] =
              True - self.columns is created
              False - Nothing happens
        I [Boolean] =
              True = a Stokes_I column will be added to the ColDefs object
              False = Stokes_I column is left out of the ColDefs object
        Q [Boolean] =
              True = a Stokes_Q column will be added to the ColDefs object
              False = Stokes_Q column is left out of the ColDefs object
        U [Boolean] =
              True = a Stokes_U column will be added to the ColDefs object
              False = Stokes_U column is left out of the ColDefs object
        V [Boolean] =
              True = a Stokes_V column will be added to the ColDefs object
              False = Stokes_V column is left out of the ColDefs object
        continuum [Boolean] =
              True = a Continuum column will be added to the ColDefs object
              False = Continuum column is left out of the ColDefs object
        """
        
        self.wl = numpy.array(self.wl)
        
        if prepareColumns:
            coldefs = []
            wave = pyfits.Column(name='Wavelength', format='D', array=self.wl)
            coldefs.append(wave)
            if I:
                flux_I = pyfits.Column(name='Stokes_I', format='D', 
                        array=numpy.array(self.flux_I))
                coldefs.append(flux_I)
            if dI:
                dflux_I = pyfits.Column(name='dStokes_I', format='D',
                        array=numpy.array(self.dflux_I))
                coldefs.append(dflux_I)
            if Q:
                flux_Q = pyfits.Column(name='Stokes_Q', format='D', 
                        array=numpy.array(self.flux_Q))
                coldefs.append(flux_Q)
            if U:
                flux_U = pyfits.Column(name='Stokes_U', format='D',
                        array=numpy.array(self.flux_U))
                coldefs.append(flux_U)
            if V:
                flux_V = pyfits.Column(name='Stokes_V', format='D',
                        array=numpy.array(self.flux_V))
                coldefs.append(flux_V)
            if continuum:
                continuum = pyfits.Column(name='Continuum', format='D',
                        array=numpy.array(self.continuum))
                coldefs.append(continuum)
            self.columns = pyfits.ColDefs(coldefs)

    def plot(self, I=True, Q=False, U=False, V=False,
                 continuum=False, ax=None, **kwargs):
        """
        Spectrum.plot(I=True, Q=False, U=False, V=False, continuum=False,
            ax=pyplot.axis, **kwargs)
            
        plot allows a simple way to plot the contents of the spectrum.
        
        I [Boolean] = signifying whether or not Stokes_I is plotted
        Q [Boolean] = signifying whether or not Stokes_Q is plotted
        U [Boolean] = signifying whether or not Stokes_U is plotted
        V [Boolean] = signifying whether or not Stokes_V is plotted
        continuum [Boolean] = signifying whether or not the continuum is plotted
        ax [matplotlib.pyplot.axis object]
        **kwargs = arguments to pass to the plot command
        """
        plotLabel = self.header.get('EXTNAME')
        if I:
            if self.dflux_I != None:
                ax.errorbar(self.wl, self.flux_I, yerr=self.dflux_I, 
                            label=plotLabel, **kwargs)
            else:
                ax.plot(self.wl, self.flux_I, label=plotLabel, **kwargs)
        if Q:
            ax.plot(self.wl, self.flux_Q, label=plotLabel, **kwargs)
        if U:
            ax.plot(self.wl, self.flux_U, label=plotLabel, **kwargs)
        if V:
            ax.plot(self.wl, self.flux_V, label=plotLabel, **kwargs)
        if continuum:
            ax.plot(self.wl, self.continuum, label=plotLabel, **kwargs)
    
    def resample(self, R=0.0, nyquist=False, observedWl=None, pad=None):
        """
        Spectrum.resample(R=0.0, nyquist=False, observedWl=None, pad=None)
        
        resample convolves the spectrum to a resolution R, and optionally
        nyquist samples it, or re-bins it to a different wavelength range,
        padding the beginning and end of the wavelength range.
        
        R [float] = Desired resolving power (Lambda/dLambda)
        nyquist [Boolean] = If true, returns a nyquist-sampled spectrum
        observedWl [numpy.array] wavelength points to which the spectrum
             will be rebinned/interpolated
        pad [None/float] If the spectrum is to be rebinned, pad contains
             the value to be stored in the flux points which are not spanned
             by a complete bin.
        """
        subsample = 16.0

        newWl = [self.wl[0]]
        while True:
            stepsize = newWl[-1]/(R*subsample)
            if newWl[-1]+stepsize > self.wl[-1]:
                break
            newWl.append(newWl[-1]+stepsize)

        if self.flux_I != None:
            I = scipy.interpolate.interpolate.interp1d(self.wl,
                    self.flux_I, bounds_error=False, fill_value=1.0)
            newI = I(newWl)
        if self.flux_V != None:
            V = scipy.interpolate.interpolate.interp1d(self.wl,
                    self.flux_V, bounds_error=False, fill_value=1.0)
            newV = V(newWl)
        if self.flux_Q != None:
            Q = scipy.interpolate.interpolate.interp1d(self.wl,
                    self.flux_Q, bounds_error=False, fill_value=1.0)
            newQ = Q(newWl)
        if self.flux_U != None:
            U = scipy.interpolate.interpolate.interp1d(self.wl,
                    self.flux_U, bounds_error=False, fill_value=1.0)
            newU = U(newWl)
        if self.continuum != None:
            continuum = scipy.interpolate.interpolate.interp1d(self.wl,
                    self.continuum, bounds_error=False, fill_value=1.0)
            newContinuum = continuum(newWl)
        const = numpy.ones(len(newWl))

        xk = numpy.array(range(int(4.0*subsample)))
        yk = numpy.exp(-(xk-(2.0*subsample))**2.0/(subsample**2.0/(4.0*numpy.log(2.0))))
    
        newWl = numpy.array(newWl[int(len(xk)/2.0):-int(len(xk)/2.0)])

        normal = scipy.signal.convolve(const, yk, mode = 'same')
        if self.flux_I != None:
            result_I = scipy.signal.convolve(newI, yk, mode ='same')/normal
            flux_I = numpy.array(result_I[int(len(xk)/2.0):-int(len(xk)/2.0)])
        else:
            flux_I = None
        if self.flux_V != None:
            result_V = scipy.signal.convolve(newV, yk, mode ='same')/normal
            flux_V = numpy.array(result_V[int(len(xk)/2.0):-int(len(xk)/2.0)])
        else:
            flux_V = None
        if self.flux_Q != None:
            result_Q = scipy.signal.convolve(newQ, yk, mode ='same')/normal
            flux_Q = numpy.array(result_Q[int(len(xk)/2.0):-int(len(xk)/2.0)])
        else:
            flux_Q = None
        if self.flux_U != None:
            result_U = scipy.signal.convolve(newU, yk, mode ='same')/normal
            flux_U = numpy.array(result_U[int(len(xk)/2.0):-int(len(xk)/2.0)])
        else:
            flux_U = None
        if self.continuum != None:
            result_Cont = scipy.signal.convolve(newContinuum, yk, mode ='same')/normal
            continuum = numpy.array(result_Cont[int(len(xk)/2.0):-int(len(xk)/2.0)])
        else:
            continuum = None

        header = self.header.copy()
        header.set('RESOLVING_POWER', R)
        processed = Spectrum(wl=newWl, I=flux_I, Q=flux_Q, U=flux_U, continuum=continuum,
                V=flux_V, header=header, spectrum_type='CONVOLVED')
                
        
        if nyquist:
            nyquistWl = []
            deltaWl = min(self.wl)/(2.0*R)
            nyquistWl.append(min(self.wl) + deltaWl)
            while True:
                deltaWl = nyquistWl[-1]/(2.0*R)
                if nyquistWl[-1]+deltaWl > self.wl[-1]:
                    break
                nyquistWl.append(nyquistWl[-1]+deltaWl)
            nyquistWl = numpy.array(nyquistWl)
            processed.bin(nyquistWl)

        if observedWl == None:
            return processed
        else:
            processed.bin(observedWl, pad=pad)
            return processed

    def rv(self, rv=0.0):
        """
        Spectrum.rv(rv=0.0)

        This routine simulates the effect of a radial velocity on the target spectrum

        rv = radial velocity in km/s

        """
        beta = rv/299792.0
        self.wl = self.wl * (1 + beta)/(1.0 - beta**2.0)**(-0.5)

    def bin(self, newWl, pad=None):
        """
        Spectrum.bin(newWl=[], pad=None)
        
        This routine simulates the binning of a synthetic spectra due to
        the discrete nature of detector pixels.
        
        newWl [numpy.array] the new wavelengths to which the spectrum
             should be binned.
        pad [None/float] pad contains the value to be stored in the flux points
             which are not spanned by a complete bin.
        """
        deltaWl = numpy.median(numpy.diff(newWl))/5.0
        if pad == None:
            interpWl = numpy.arange(self.wl[0], self.wl[-1], deltaWl)
        else:
            interpWl = numpy.arange(newWl[0], newWl[-1], deltaWl)
        newWave = []

        if self.flux_I != None:
            I = scipy.interpolate.splrep(self.wl, self.flux_I)
            I_interp = scipy.interpolate.splev(interpWl, I, ext=1)
            if pad != None:
                I_interp[I_interp==0]=1.0
            newSpec_I = []
        if self.flux_Q != None:
            Q = scipy.interpolate.splrep(self.wl, self.flux_Q)
            Q_interp = scipy.interpolate.splev(interpWl, Q, ext=1)
            newSpec_Q = numpy.zeros(len(newWl))
        if self.flux_U != None:
            U = scipy.interpolate.splrep(self.wl, self.flux_U)
            U_interp = scipy.interpolate.splev(interpWl, U, ext=1)
            newSpec_U = numpy.zeros(len(newWl))
        if self.flux_V != None:
            V = scipy.interpolate.splrep(self.wl, self.flux_V)
            V_interp = scipy.interpolate.splev(interpWl, V, ext=1)
            if pad != None:
                V_interp[V_interp==0]=0.0
            newSpec_V = []
        if self.continuum != None:
            continuum = scipy.interpolate.splrep(self.wl, self.continuum)
            cont_interp = scipy.interpolate.splev(interpWl, continuum, ext=1)
            newSpec_continuum = numpy.zeros(len(newWl))
        for i in range(len(newWl)):
            if i==0:
                lowerBound = newWl[0]-deltaWl*5.0
            else:
                lowerBound = (newWl[i-1]+newWl[i])/2.0
            if i==len(newWl)-1:
                upperBound = newWl[-1]+deltaWl*5.0
            else:
                upperBound = (newWl[i]+newWl[i+1])/2.0
            inBin = scipy.where( (interpWl > lowerBound) & (
                interpWl <= upperBound))[0]
            if (len(inBin) > 1):
                newWave.append(newWl[i])
                denom = interpWl[inBin][-1] - interpWl[inBin][0]
                if self.flux_I != None:
                    num=scipy.integrate.simps(I_interp[inBin], 
                            x=interpWl[inBin])
                    newSpec_I.append(num/denom)
                if self.flux_Q != None:
                    num=scipy.integrate.simps(Q_interp[inBin], 
                            x=interpWl[inBin])
                    newSpec_Q[i] = num/denom
                if self.flux_U != None:
                    num=scipy.integrate.simps(U_interp[inBin], 
                            x=interpWl[inBin])
                    newSpec_U[i] = num/denom
                if self.flux_V != None:
                    num=scipy.integrate.simps(V_interp[inBin], 
                            x=interpWl[inBin])
                    newSpec_V.append(num/denom)
                if self.continuum != None:
                    num=scipy.integrate.simps(cont_interp[inBin], 
                            x=interpWl[inBin])
                    newSpec_continuum[i] = num/denom

        self.wl = numpy.array(newWave)
        if self.flux_I != None:
            self.flux_I = numpy.array(newSpec_I)
        if self.flux_Q != None:
            self.flux_Q = numpy.array(newSpec_Q)
        if self.flux_U != None:
            self.flux_U = numpy.array(newSpec_U)
        if self.flux_V != None:
            self.flux_V = numpy.array(newSpec_V)
        if self.continuum != None:
            self.continuum = numpy.array(newSpec_continuum)

    def rotate(self, angle=0.0, wlPoint = None):
        """
        Spectrum.rotate(angle=0.0, wlPoint = None)
        
        This routine rotates the spectra, to mimic errors in the continuum
             determination.
        angle [rad] = arctan(rise/run).

        Units are continuum/angstrom
        """
        
        I = None
        Q = None
        U = None
        V = None
        continuum = None

        if wlPoint == None:
            wlPoint = (self.wl[0]+self.wl[-1])/2.0
        
        offset = numpy.tan(angle)*(self.wl-wlPoint)
        if (self.flux_I != None):
            I = self.flux_I + offset
        if (self.flux_Q != None):
            Q = self.flux_Q + offset
        if (self.flux_U != None):
            U = self.flux_U + offset
        if (self.flux_V != None):
            V = self.flux_V + offset
        if (self.continuum != None):
            continuum = self.continuum + offset

        return Spectrum(wl=self.wl, I=I, Q=Q, U=U, V=V, 
                continuum=continuum, header=self.header,
                spectrum_type="ROTATED")


    def __sub__(self, other):
        '''
        Spectrum.__sub__(other)
		
        __sub__ overloads the subtraction operator.  It subtracts one spectrum 
        from the other
        
        subtracted = Spectrum - other
        '''
        overlap_start = numpy.max([numpy.min(self.wl), numpy.min(other.wl)])
        overlap_stop = numpy.min([numpy.max(self.wl), numpy.max(other.wl)])
        overlap_self = scipy.where((self.wl >= overlap_start) & (self.wl <= overlap_stop))
        overlap_other = scipy.where((other.wl >= overlap_start) & (other.wl <= overlap_stop))
        
        I = None
        Q = None
        U = None
        V = None
        continuum = None
        
        if (self.flux_I != None) & (other.flux_I != None):
            I = self.flux_I[overlap_self] - other.flux_I[overlap_other]
        if (self.flux_Q != None) & (other.flux_Q != None):
            Q = self.flux_Q[overlap_self] - other.flux_Q[overlap_other]
        if (self.flux_U != None) & (other.flux_U != None):
            U = self.flux_U[overlap_self] - other.flux_U[overlap_other]
        if (self.flux_V != None) & (other.flux_V != None):
            V = self.flux_V[overlap_self] - other.flux_V[overlap_other]
        if (self.continuum != None) & (other.continuum != None):
            continuum = self.continuum[overlap_self] - other.continuum[overlap_other]
        
        return Spectrum(wl=self.wl[overlap_self], I=I, Q=Q, U=U, V=V, 
                continuum=continuum, header=self.header,
                spectrum_type="DIFFERENCE")

    def __div__(self, factor):
        """
        Spectrum.__div__(factor)
        
        __div__ overloads the division operator.  If factor is a scalar, div
                  returns a spectrum object divided by the scalar factor.  If
                  factor is instead another Spectrum object, div returns a
                  Spectrum object of one spectrum divided by the other.
        
        divided = Spectrum/factor
        """

        if isinstance(factor, float):
            I = None
            Q = None
            U = None
            V = None
            continuum = None

            if (self.flux_I != None):
                I = self.flux_I/factor
            if (self.flux_Q != None):
                Q = self.flux_Q/factor
            if (self.flux_U != None):
                U = self.flux_U/factor
            if (self.flux_V != None):
                V = self.flux_V/factor
            if (self.continuum != None):
                continuum = self.continuum/factor

            return Spectrum(wl=self.wl, I=I, Q=Q, U=U, V=V, continuum=continuum,
                    header=self.header, spectrum_type="SCALED")

        elif isinstance(factor, Spectrum):
            overlap_start = numpy.max([numpy.min(self.wl), numpy.min(factor.wl)])
            overlap_stop = numpy.min([numpy.max(self.wl), numpy.max(factor.wl)])
            overlap = scipy.where((self.wl >= overlap_start) & (self.wl <= overlap_stop))

            if (self.flux_I != None) & (factor.flux_I != None):
                I = scipy.interpolate.splrep(factor.wl, factor.flux_I)
                retval_I = numpy.zeros(len(self.wl))
                retval_I[overlap] = self.flux_I[overlap] - scipy.interpolate.splev(self.wl[overlap],I)
            else:
                retval_I = None
            if (self.flux_Q != None) & (factor.flux_Q != None):
                Q = scipy.interpolate.splrep(factor.wl, factor.flux_Q)
                retval_Q = numpy.zeros(len(self.wl))
                retval_Q[overlap] = self.flux_Q[overlap] - scipy.interpolate.splev(self.wl[overlap],Q)
            else:
                retval_Q = None
            if (self.flux_U != None) & (factor.flux_U != None):
                U = scipy.interpolate.splrep(factor.wl, factor.flux_U)
                retval_U = numpy.zeros(len(self.wl))
                retval_U[overlap] = self.flux_U[overlap] - scipy.interpolate.splev(self.wl[overlap],U)
            else:
                retval_U = None
            if (self.flux_V != None) & (factor.flux_V != None):
                V = scipy.interpolate.splrep(factor.wl, factor.flux_V)
                retval_V = numpy.zeros(len(self.wl))
                retval_V[overlap] = self.flux_V[overlap] - scipy.interpolate.splev(self.wl[overlap],V)
            else:
                retval_V = None
            if (self.continuum != None) & (factor.continuum != None):
                continuum = scipy.interpolate.splrep(factor.wl, factor.continuum)
                retval_continuum = numpy.zeros(len(self.wl))
                retval_continuum[overlap] = self.continuum[overlap] - scipy.interpolate.splev(self.wl[overlap],continuum)
            else:
                retval_continuum = None
            return Spectrum(wl=self.wl, I=retval_I, Q=retval_Q, U=retval_U, V=retval_V,
                    continuum=retval_continuum, header=self.header, 
                    spectrum_type="DIVIDED")

    def diff_spectra(self, other, pad=False):
        '''
        Spectrum.diff_spectrum(other, pad=False)
        
        diff_spectrum computes the difference between this spectrum and another
        
        other = Spectrum object containing the comparison spectrum
        pad = Boolean determining whether or not the difference spectrum should
            be zero-padded
        '''
        if self.R != other.R:
            errmsg = "The resolutions of the two spectra are not compatible!"
            raise SpectrumError(4, errmsg)

        overlap_start = numpy.max([numpy.min(self.wl), numpy.min(other.wl)])
        overlap_stop = numpy.min([numpy.max(self.wl), numpy.max(other.wl)])
        overlap = scipy.where((self.wl >= overlap_start) & (self.wl <= overlap_stop))

        """
        if (self.flux_I != None) & (other.flux_I != None):
            I = scipy.interpolate.splrep(other.wl, other.flux_I)
        if (self.flux_Q != None) & (other.flux_Q != None):
            Q = scipy.interpolate.splrep(other.wl, other.flux_Q)
        if (self.flux_U != None) & (other.flux_U != None):
            U = scipy.interpolate.splrep(other.wl, other.flux_U)
        if (self.flux_V != None) & (other.flux_V != None):
            V = scipy.interpolate.splrep(other.wl, other.flux_V)
        if (self.continuum != None) & (other.continuum != None):
            continuum = scipy.interpolate.splrep(other.wl, other.continuum)
        """
        if pad:
            if (self.flux_I != None) & (other.flux_I != None):
                I = scipy.interpolate.splrep(other.wl, other.flux_I)
                retval_I = numpy.zeros(len(self.wl))
                retval_I[overlap] = self.flux_I[overlap] - scipy.interpolate.splev(self.wl[overlap],I)
            else:
                retval_I = None
            if (self.flux_Q != None) & (other.flux_Q != None):
                Q = scipy.interpolate.splrep(other.wl, other.flux_Q)
                retval_Q = numpy.zeros(len(self.wl))
                retval_Q[overlap] = self.flux_Q[overlap] - scipy.interpolate.splev(self.wl[overlap],Q)
            else:
                retval_Q = None
            if (self.flux_U != None) & (other.flux_U != None):
                U = scipy.interpolate.splrep(other.wl, other.flux_U)
                retval_U = numpy.zeros(len(self.wl))
                retval_U[overlap] = self.flux_U[overlap] - scipy.interpolate.splev(self.wl[overlap],U)
            else:
                retval_U = None
            if (self.flux_V != None) & (other.flux_V != None):
                V = scipy.interpolate.splrep(other.wl, other.flux_V)
                retval_V = numpy.zeros(len(self.wl))
                retval_V[overlap] = self.flux_V[overlap] - scipy.interpolate.splev(self.wl[overlap],V)
            else:
                retval_V = None
            if (self.continuum != None) & (other.continuum != None):
                continuum = scipy.interpolate.splrep(other.wl, other.continuum)
                retval_continuum = numpy.zeros(len(self.wl))
                retval_continuum[overlap] = self.continuum[overlap] - scipy.interpolate.splev(self.wl[overlap],continuum)
            else:
                retval_continuum = None
            return Spectrum(wl=self.wl, I=retval_I, Q=retval_Q, U=retval_U,
                   V=retval_V, continuum=retval_continuum, header=self.header,
                   spectrum_type='DIFFERENCE')
        else:
            """
            To be consistent, this case should return a Spectrum object.
            """
            return numpy.array(x1)[overlap], numpy.array(y1)[overlap] - scipy.interpolate.splev(x1[overlap],y)

    def blend(self, other, fraction):
        """
        blended = Spectrum.blend(other, fraction)
        
        Spectrum.blend returns a linear blend of the current spectrum with the other
              Spectrum object weighted by the scalar fraction
              
        other : [Spectrum] - the other spectrum
        fraction : [float] - the ratio of blending, obeying the limits:
                    0 - all other spectrum
                    0.5 - equal blend
                    1 - all this spectrum

        The function returns a Spectrum object containing the blended spectrum
        """
        overlap_start = numpy.max([numpy.min(self.wl), numpy.min(other.wl)])
        overlap_stop = numpy.min([numpy.max(self.wl), numpy.max(other.wl)])
        overlap = scipy.where((self.wl >= overlap_start) & (self.wl <= overlap_stop))

        newWl = self.wl[overlap]
        newI = None
        newQ = None
        newU = None
        newV = None
        newCont = None

        if (self.flux_I != None) & (other.flux_I != None):
            I = scipy.interpolate.splrep(other.wl, other.flux_I)
            newI = self.flux_I[overlap]*fraction + scipy.interpolate.splev(newWl, I)*(1.0-fraction)
        if (self.flux_Q != None) & (other.flux_Q != None):
            Q = scipy.interpolate.splrep(other.wl, other.flux_Q)
            newQ = self.flux_I[overlap]*fraction + scipy.interpolate.splev(newWl, Q)*(1.0-fraction)
        if (self.flux_U != None) & (other.flux_U != None):
            U = scipy.interpolate.splrep(other.wl, other.flux_U)
            newU = self.flux_I[overlap]*fraction + scipy.interpolate.splev(newWl, U)*(1.0-fraction)
        if (self.flux_V != None) & (other.flux_V != None):
            V = scipy.interpolate.splrep(other.wl, other.flux_V)
            newV = self.flux_V[overlap]*fraction + scipy.interpolate.splev(newWl, V)*(1.0-fraction)
        if (self.continuum != None) & (other.continuum != None):
            continuum = scipy.interpolate.splrep(other.wl, other.flux_V)
            newCont = self.flux_I[overlap]*fraction + scipy.interpolate.splev(newWl, I)*(1.0-fraction)

        return Spectrum(wl=newWl, I=newI, Q=newQ, U=newU, V=newV, continuum=newCont, header=self.header,
                        spectrum_type="BLENDED")
    
    def calc_EW(self, wlStart, wlStop, findContinuum=False):
        """
        EW = Spectrum.calc_EW(wlStart, wlStop, findContinuum=False)
        
        calc_EW calculates the equivalent width of the Spectrum object between the 
             given start and stop wavelengths.
             
        wlStart [float] = start of the EW interval.  Must be same units as the
             Spectrum.wl array
        wlStop [float] = stop of the EW interval.  Must be same units as the
             Spectrum.wl array
        findContinuum [Boolean] = Whether or not to attempt to automatically find
             the continuum.  Should probably only be used if Spectrum.flux_I is 
             not normalized
        """
        if (wlStart > self.wl[-1]) or (wlStop < self.wl[0]):
            raise SpectrumError(2, 'Wavelength Regions do not overlap!')

        bm = scipy.where( (self.wl > wlStart) & (self.wl < wlStop) )[0]
        cont = numpy.ones(len(bm))
        if findContinuum:
            cont *= numpy.median(self.flux_I[bm])
            print "%.4f - continuum level" % numpy.median(self.flux_I[bm])
        num = scipy.integrate.simps(self.flux_I[bm], self.wl[bm])
        denom = scipy.integrate.simps(cont, self.wl[bm])
        return (denom-num)

    def mergeSpectra(self, second=None):
        """
        merged = Spectrum.mergeSpectra(second=None)
        
        Spectrum.mergeSpectra merges the spectrum with another spectrum which
            covers a different spectral region.
            
        second [Spectrum] = Spectrum object to be merged
        """
        if second == None:
            return self

        x1 = self.wl
        x2 = second.wl


        overlap_start = numpy.max([numpy.min(x1), numpy.min(x2)])
        overlap_stop = numpy.min([numpy.max(x1), numpy.max(x2)])
        overlap = scipy.where((x1 >= overlap_start) & (x1 <= overlap_stop))
        if (len(overlap[0]) > 1):
            unique1 = scipy.where((x1 < overlap_start) | (x1 > overlap_stop))
            unique2 = scipy.where((x2 < overlap_start) | (x2 > overlap_stop))

            new_x = numpy.append(x1, x2[unique2])
    
            if (self.flux_I != None) & (second.flux_I != None):
                I1 = self.flux_I
                I2 = second.flux_I
                I1[numpy.isnan(I1)] = 0.0
                I2[numpy.isnan(I2)] = 0.0
                I = scipy.interpolate.splrep(x2, I2)
                Iinterp = scipy.interpolate.splev(x1[overlap], I)
                if (self.dflux_I != None) & (second.dflux_I != None):
                    dI1 = self.dflux_I
                    dI2 = second.dflux_I
                    dI1[numpy.isnan(dI1)] = 1000000.0
                    dI2[numpy.isnan(dI2)] = 1000000.0
                    dI1[I1==0.0] = 100000000.0
                    dI2[I2==0.0] = 100000000.0
                    dI = scipy.interpolate.splrep(x2, dI2)
                    dIinterp = scipy.interpolate.splev(x1[overlap], dI)
                    """
                    import matplotlib.pyplot as pyplot
                    fig = pyplot.figure(1)
                    fig.clear()
                    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
                    #print dIinterp
                    print overlap_start, overlap_stop
                    #raw_input()
                    #"""
                    mergedI = (I1[overlap]/dI1[overlap] + Iinterp/dIinterp)/(1.0/dI1[overlap]+1.0/dIinterp)
                    merged_dI = 1.0/(1.0/dI1[overlap] + 1.0/dIinterp)
                    new_dI = numpy.append(numpy.append(dI1[unique1], merged_dI), dI2[unique2])
                    """
                    ax1.errorbar(x1[overlap], I1[overlap], yerr=dI1[overlap], color = 'r')
                    ax1.errorbar(x2, I2, yerr=dI2, color = 'b')
                    ax1.errorbar(x1[overlap], mergedI, yerr=merged_dI, color = 'g')
                    fig.show()
                    raw_input()
                    #"""
                else:
                    mergedI = (I1[overlap] + Iinterp)/2.0
                    new_dI = None
                new_I = numpy.append(numpy.append(I1[unique1], mergedI), I2[unique2])
            else:
                new_dI = None
                new_I = None
            if (self.flux_Q != None) & (second.flux_Q != None):
                Q1 = self.flux_Q
                Q2 = second.flux_Q
                Q1[numpy.isnan(Q1)] = 0.0
                Q2[numpy.isnan(Q2)] = 0.0
                Q = scipy.interpolate.splrep(x2, Q2)
                Qinterp = scipy.interpolate.splev(x1[overlap], Q)
                mergedQ = (Q1[overlap] + Qinterp)/2.0
                new_Q = numpy.append(numpy.append(Q1[unique1], mergedQ), Q2[unique2])
            else:
                new_Q = None
            if (self.flux_U != None) & (second.flux_U != None):
                U1 = self.flux_U
                U2 = second.flux_U
                U1[numpy.isnan(U1)] = 0.0
                U2[numpy.isnan(U2)] = 0.0
                U = scipy.interpolate.splrep(x2, U2)
                Uinterp = scipy.interpolate.splev(x1[overlap], U)
                mergedU = (U1[overlap] + Uinterp)/2.0
                new_U = numpy.append(numpy.append(U1[unique1], mergedU), U2[unique2])
            else:
                new_U = None
            if (self.flux_V != None) & (second.flux_V != None):
                V1 = self.flux_V
                V2 = second.flux_V
                V1[numpy.isnan(V1)] = 0.0
                V2[numpy.isnan(V2)] = 0.0
                V = scipy.interpolate.splrep(x2, V2)
                Vinterp = scipy.interpolate.splev(x1[overlap], V)
                mergedV = (V1[overlap] + Vinterp)/2.0
                new_V = numpy.append(numpy.append(V1[unique1], mergedV), V2[unique2])
            else:
                new_V = None
            if (self.continuum != None) & (second.continuum != None):
                C1 = self.continuum
                C2 = second.continuum
                C1[numpy.isnan(C1)] = 0.0
                C2[numpy.isnan(C2)] = 0.0
                C = scipy.interpolate.splrep(x2, C2)
                Cinterp = scipy.interpolate.splev(x1[overlap], C)
                mergedC = (C1[overlap] + Cinterp)/2.0
                new_C = numpy.append(numpy.append(C1[unique1], mergedC), C2[unique2])
            else:
                new_C = None
        else:
            new_x = numpy.append(x1, x2)
            if (self.flux_I != None) & (second.flux_I != None):
                I1 = self.flux_I
                I2 = second.flux_I
                I1[numpy.isnan(I1)] = 0.0
                I2[numpy.isnan(I2)] = 0.0
                new_I = numpy.append(I1, I2)
                if (self.dflux_I != None) & (second.dflux_I != None):
                    dI1 = self.dflux_I
                    dI2 = second.dflux_I
                    dI1[numpy.isnan(dI1)] = 0.0
                    dI2[numpy.isnan(dI2)] = 0.0
                    new_dI = numpy.append(dI1, dI2)
                else:
                    new_dI = None
            else:
                new_I = None
                new_dI = None
            if (self.flux_Q != None) & (second.flux_Q != None):
                Q1 = self.flux_Q
                Q2 = second.flux_Q
                Q1[numpy.isnan(Q1)] = 0.0
                Q2[numpy.isnan(Q2)] = 0.0
                new_Q = numpy.append(Q1, Q2)
            else:
                new_Q = None
            if (self.flux_U != None) & (second.flux_U != None):
                U1 = self.flux_U
                U2 = second.flux_U
                U1[numpy.isnan(U1)] = 0.0
                U2[numpy.isnan(U2)] = 0.0
                new_U = numpy.append(U1, U2)
            else:
                new_U = None
            if (self.flux_V != None) & (second.flux_V != None):
                V1 = self.flux_V
                V2 = second.flux_V
                V1[numpy.isnan(V1)] = 0.0
                V2[numpy.isnan(V2)] = 0.0
                new_V = numpy.append(V1, V2)
            else:
                new_V = None
            if (self.continuum != None) & (second.continuum != None):
                C1 = self.continuum
                C2 = second.continuum
                C1[numpy.isnan(C1)] = 0.0
                C2[numpy.isnan(C2)] = 0.0
                new_C = numpy.append(C1, C2)
            else:
                new_C = None
            
        header = self.header.copy()
        header.set("WLSTART", numpy.min([header.get("WLSTART"), 
            second.header.get("WLSTART")]))
        header.set("WLSTOP", numpy.max([header.get("WLSTOP"), 
            second.header.get("WLSTOP")]))


        retval = Spectrum(wl=new_x, I = new_I, dI = new_dI, Q = new_Q, U = new_U, V = new_V, 
                continuum = new_C, spectrum_type='MERGED', header=header)

        return retval

class ObservedSpectrum ( object ):
    def __init__(self, observed=None):
        self.observed = observed             # Should be a Spectrum object

    def yank(self, **kwargs):
        return self.observed

class Integrator( object ):
    def __init__(self, parent=None, deltav = 0.1, limb_darkening=None):
        """
        Integrator(parent=None, deltav=0.1, limb_darkening=None)
            
        An Integrator object handles the interpolating, integrating, and
            convolving of raw data produced by MoogStokes
        
        Input:
            parent [Moog960.SyntheticPhrase] = reference to the Integrator's parent
            deltav [float] = velocity resolution used for interpolating (km/s)
            limb_darkening [numpy.array(float)] = custom limb darkening coefficients for
                the disk integration routines.  There should be as many coefficients as 
                there are unique emergent spectral elements.
        Contents:
        
        Integrator.parent = reference to Integrator's parent.  Necessary to find the
            location of the raw/processed data
        Integrator.deltav = velocity resolution used for interpolating (km/s)
        Integrator.interpolatedData = reference to interpolated data
        Integrator.integratedData = reference to integrated data
        Integrator.convolvedData = reference to convolved data
        Integrator.limb_darkening = custom limb darkening coefficients for disk integration 
            routines	
        """
        self.parent = parent
        self.deltav = deltav
        self.interpolated = parent.interpolatedData  # interpolated to uniform wl grid
        for interp in self.interpolated:
            interp.loadData()
        self.integrated = parent.integratedData      # vsin i - needs interpolated
        for integ in self.integrated:
            integ.loadData()
        self.convolved = parent.convolvedData        # R - needs integrated
        for convol in self.convolved:
            if convol.wl == None:
                convol.loadData()
        self.limb_darkening = limb_darkening

class TennisBall( Integrator ):
    """
    TennisBall - inherits from Integrator
    
    TennisBall is used to handle data processing from the normal, non-magnetic
    Moog.  Tennis balls are all one color (I know, I know, except for the seam)
    and the scalar version of Moog produces disk-averaged spectra.  There is
    no need to actually process the output of Moog, so the diskInt, resample
    """
    def loadData(self):
        """
        TennisBall.loadData()
        
        placeholder function.  Since no data processing is necessary,
        nothing is done.
        """
        return

    def diskInt(self):
        """
        TennisBall.diskInt()
        
        Copies the first entry in the parent.rawData list to the list of 
        integrated spectra.
        
        Since the TennisBall object is used with Moog-produced disk-averaged
        spectra, no processing is performed.
        """
        self.integrated.append(self.parent.rawData[0])

    def resample(self, R=0, observedWl=None):
        """
        TennisBall.resample(R=0, observedWl=None)
        
        Resamples the first entry in the list of integrated spectra and 
        saves the result in the list of convolved spectra.
        """
        self.convolved.append(self.integrated[0].resample(R=R, observedWl=observedWl))

    def yank(self, vsini=0.0, R=0.0, observedWl = None, keySignature="CONVOLVED"):
        """
        TennisBall.yank(vsini=0.0, R=0.0, observedWl = None, 
               keySignature= "RAW", "INTEGRATED", "CONVOLVED")
        """
        if keySignature == "CONVOLVED":
            return self.convolved[0]
        if keySignature == "RAW":
            return self.parent.rawData[0]
        if keySignature == "INTEGRATED":
            return self.integrated[0]

class BeachBall( Integrator ):
    """
    A BeachBall object is an extension of the Integrator class and is used to
    handle MoogStokes output generated with the MoogStokes parameter diskflag=1
    
    The case of the BeachBall disk integration algorithm divides the stellar surface into
    N 
    """
    def loadData(self):
        """
        BeachBall.loadData()
        
        Loads raw data from the parent Moog960.SyntheticPhrase.  Then, 
        interpolates each emergent spectra to the wavelength spacing denoted
        by the velocity resolution (BeachBall.deltav), weighting each
        slice by a limb darkening coefficient.  The interpolated and weighted
        emergent fluxes are then stored in the interpolated list
        """
        c = 3e5              #km/s
        phi = []
        mu = []
        cell = []

        newWl = []
        newWl.append(self.parent.rawData[0].wl[0])
        while newWl[-1] < self.parent.rawData[0].wl[-1]:
            dLambda = newWl[-1]*self.deltav/c
            newWl.append(newWl[-1]+dLambda)
        newWl = numpy.array(newWl[0:-1])

        wave = numpy.mean(newWl)
        if ((1.0/(wave/10000.0)) < 2.4):
            self.alpha = -0.023 + 0.292/(wave/10000.0)
        else:
            self.alpha = -0.507 + 0.441/(wave/10000.0)

        limb_darkening = []
        for raw in self.parent.rawData:
            phi.append(raw.header.get('PHI_ANGLE'))
            mu.append(raw.header.get('MU'))
            limb_darkening.append(1.0-(1.0-mu[-1]**(self.alpha)))
            cell.append(raw.header.get('CELL'))
            fI = scipy.interpolate.UnivariateSpline(raw.wl, 
                    numpy.array(raw.flux_I)/numpy.array(raw.continuum), s=0)
            fV = scipy.interpolate.UnivariateSpline(raw.wl, 
                    numpy.array(raw.flux_V)/numpy.array(raw.continuum), s=0)
            self.interpolated.append(Spectrum(wl=newWl, I = fI(newWl)*limb_darkening[-1],
                V = fV(newWl)*limb_darkening[-1], continuum = numpy.ones(len(newWl))*limb_darkening[-1],
                header = raw.header.copy(), spectrum_type='INTERPOLATED'))
                
            """
            Probably should do something with the labels here...
            """

        self.limb_darkening = numpy.array(limb_darkening)
        self.phi = numpy.array(phi)
        self.mu = numpy.array(mu)
        self.cell = numpy.array(cell)
        self.ncells = len(self.cell)


    def diskInt(self, vsini=0.0):
        """
        BeachBall.diskInt(vsini=0.0)
        
        Input:
            vsini [float] = rotational velocity * sin (inclination) (in km/s)
        
        Returns:
            True:       if the requested VSINI disk integrated spectrum
                        did not already exist and had to be created
            False:      if the requested VSINI disk integrated spectrum was
                        already existing, and did NOT have to be created
        """
        if (self.parent.rawData[0].wl == None):
            for raw in self.parent.rawData:
                raw.loadData()
            self.loadData()
        for integrated in self.integrated:
            if integrated.header.get('VSINI') == vsini:
                return False
        I, V = self.rtint(vsini_in=vsini)
        header = self.interpolated[0].header.copy()
        header.set('VSINI', vsini)
        header.remove('PHI_ANGLE')
        header.remove('MU')
        header.remove('CELL')
        for interp in self.interpolated[1:]:
            header.add_history(interp.header.get('SPECTRUM_TYPE')+' - '+interp.header.get('SPECTRUM_ID'))

        self.integrated.append(Spectrum(wl=self.interpolated[0].wl, I=I, V=V, header=header,
             spectrum_type='DISK INTEGRATED'))
        return True
        
        """
        Probably should do something with the labels here
        """

    #"""
    def findVsini(self, vsini):
        for integrated in self.integrated:
            if numpy.abs(integrated.header.get('VSINI') - vsini) < 0.01:
                return integrated

        raise SpectrumError(1, "Integrated Spectrum with vsini=%.2f NOT FOUND!!!" % 
                (vsini))
        #self.diskInt(vsini=vsini)
        #return self.integrated[-1]
    #"""

    def resample(self, vsini=0.0, R=0, observedWl=None):
        """
        BeachBall.resample(vsini=0.0, R=0, observedWl=None)
        
        Resample resamples the current spectrum to a desired resolving power
        
        Returns:
            retval : List of strings containing types of spectra created
                     "INTEGRATED" - If disk-integrated spectra with desired VSINI
                                does not exist, it must be created
                     "CONVOLVED" - If 
                     
        This should be done with the Labels.
        """
        retval = []
        for convol in self.convolved:
            if ((numpy.abs(convol.header.get('VSINI') - vsini) < 0.01)
             and (numpy.abs(convol.header.get('RESOLVING_POWER') - R) < 0.1)):
                return retval

        try:
            integrated = self.findVsini(vsini)
        except SpectrumError:
            self.diskInt(vsini=vsini)
            retval.append("INTEGRATED")
            integrated = self.findVsini(vsini)
        if R > 0:
            self.convolved.append(integrated.resample(R, observedWl=observedWl))
            retval.append("CONVOLVED")
            return retval
        else:
            raise SpectrumError(3, "Resolving Power must be greater than 0!")

    def yank(self, vsini=0.0, R=0.0, observedWl = None, keySignature="CONVOLVED"):
        if keySignature=="INTEGRATED":
            return self.findVsini(vsini)

        if keySignature=="CONVOLVED":
            for convol in self.convolved:
                if ((numpy.abs(convol.header.get('VSINI') - vsini) < 0.01) and
                    (numpy.abs(convol.header.get('RESOLVING_POWER') - R) < 0.1)):
                    if observedWl!= None:
                        return rebin(convol, observedWl)
                    else:
                        return convol

        raise SpectrumError(1, "Spectrum with vsini=%.2f and R=%.1f NOT FOUND!!!" % 
                (vsini, R))


    def rtint(self, vsini_in=0.0, vrt_in=0, **kwargs):
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
        mu = self.mu

        if "OSAMP" in kwargs:
            os = max(round(kwargs["OSAMP"]), 1)
        else:
            os = 1

        #Convert input MU to proj. radii, R of annuli for star of unit radius
        #(which is just sine rather than cosine of the angle between the outward
        #normal and the LOS)
        rmu = numpy.sqrt(1.0-self.mu**2)

        #Sort the proj. radii and corresponding intensity spectra into ascending
        #order (i.e. from disk center to limb), which is equivalent to sorting
        #MU in decending order
        order = numpy.argsort(rmu)
        rmu = rmu[order]
        nmu = len(self.mu)
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

        npts = len(self.interpolated[0].flux_I)
        xpix = numpy.arange(npts)
        nfine = os*npts
        xfine = 0.5/os * 2.0*numpy.arange(nfine)-os+1

        #Loop through annuli, constructing and convolving with rotation kernels.
        dummy = 0
        Ifine = numpy.zeros(nfine)
        Vfine = numpy.zeros(nfine)
        cfine = numpy.zeros(nfine)
        fluxI = numpy.zeros(nfine)
        fluxV = numpy.zeros(nfine)
        continuum = numpy.zeros(nfine)
        for m, spectrum, w, i in zip(mu, self.interpolated, wt, range(nmu)):
            I = spectrum.flux_I
            V = spectrum.flux_V
            c = spectrum.continuum
            #use cubic spline routine to make an oversampled version of the
            #intensity profile for the current annulus.
            if os== 1:
                Ifine = I.copy()
                Vfine = V.copy()
                cfine = c.copy()
            else:
                Ispl = scipy.interpolate.splrep(xpix, I)
                Vspl = scipy.interpolate.splrep(xpix, V)
                cspl = scipy.interpolate.splref(xpix, c)
                Ifine = scipy.interpolate.splev(Ispl, xfine)
                Vfine = scipy.interpolate.splev(Vspl, xfine)
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
                dv = self.deltav/os
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
                #print("dv = %.2f" % dv)
                #print("vsini = %.2f" % vsini)
                #print abs(v)
                #print vsini*r1
                #print vsini*r2
                #print("len(j1) = %d len(j2) = %d" % (len(j1[0]), len(j2[0])))
                #print("Rkern = ")
                #print rkern
                rkern = rkern / rkern.sum()   # normalize kernel


        # Convolve the intensity profile with the rotational velocity kernel for
        # this annulus.  Pad end of each profile with as many points as are in
        # the convolution kernel, reducing Fourier ringing.  The convolution 
        # may also be done with a routine called "externally" which efficiently
        # shifts and adds.
                if nrk > 3:
                    Ifine = scipy.convolve(Ifine, rkern, mode='same')
                    Vfine = scipy.convolve(Vfine, rkern, mode='same')
                    cfine = scipy.convolve(cfine, rkern, mode='same')

        # Calc projected simga for radial and tangential velocity distributions.
            sigma = os*vrt/numpy.sqrt(2.0) /self.deltav
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
            Ifine = scipy.convolve(Ifine, mkern, mode='same')
            Vfine = scipy.convolve(Vfine, mkern, mode='same')
            cfine = scipy.convolve(cfine, mkern, mode='same')

    # Add contribution from current annulus to the running total
            fluxI += w*Ifine
            fluxV += w*Vfine
            continuum += w*cfine

        return fluxI/continuum, fluxV/continuum




#class Diskoball( Integrator ):


"""
class MoogStokesSpectrum( Spectrum ):
    def __init__(self, name='', memory=False, **kwargs):
        self.memory = memory
        self.name = name
        if self.memory:
            self.parent = kwargs["PARENT"]
            self.deltav = self.parent.deltav
            self.vsini = self.parent.vsini
        else:
            self.deltav = kwargs["DELTAV"]
            self.vsini = kwargs["VSINI"]

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
            I = numpy.array(self.parent.flux_I)/numpy.array(self.parent.continuum)
            V = numpy.array(self.parent.flux_V)/numpy.array(self.parent.continuum)
            self.continuum = numpy.array(self.parent.continuum)
            self.wl = numpy.array(self.parent.wave)
        else:
            if "wlRange" in kwargs.keys():
                self.wlStart = kwargs["wlRange"][0]
                self.wlStop = kwargs["wlRange"][1]

    #def saveRaw(self):
#"""


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
        fileNames = ['U_Landolt.dat', 'B_Bessell.dat', 'V_Bessell.dat', 'cousins_Rband.dat', 
                    'cousins_Iband.dat', 'J_2MASS.dat', 'H_2MASS.dat', 'K_2MASS.dat']
        fnu_zero = [1829, 4144, 3544, 2950, 2280.0, 1594.0, 1024.0, 666.7 ]
        flam_zero = [4.0274905e-09, 6.3170333e-09, 3.6186341e-09, 2.1651655e-9, 1.1326593e-09, 
                  3.129e-10, 1.133e-10, 4.283e-11] #erg/s/cm^2/Angstrom
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
