from flipper import *
import scipy.ndimage

def optimalFilter(mp, noiseSpectrum = None, signalTransform = None, fwhm=1.4, extraFilt = None, 
        forceGaussian = False):
    """
    @param mp liteMap to filter
    @param noiseSpectrum an array containing the noise power spectrum of mp
    @param signalTransform an array containing the fourier transform of the signal for matched filtering
    @param fwhm if signalTransform is not specified, then filter with a gaussian of this width
    @param extraFilt another filter array to apply (form = liteMap)
    """

    """ get power spectra of cmb and noise maps """
    pNoise = liteMap.fftTools.powerFromLiteMap(mp)
    if noiseSpectrum != None:
        pNoise.powerMap = noiseSpectrum
    
    """ calculate filter """
    filt = 1 / pNoise.powerMap

    if noiseSpectrum == None: #if we are not using a nice model, the filter need conditioning...
        #eliminate outliers from filter
        med = numpy.median(filt)
        filt[numpy.where(filt> 10*med)] = med
        # convolve filter with gaussian to smooth it
        kernel_width = (5,5)  # width of gaussian kernal in sigma
        filt = scipy.ndimage.gaussian_filter(filt, kernel_width)

    if extraFilt != None:
        filt *= extraFilt
    
    ft = liteMap.fftFromLiteMap(mp)
    u = ft.kMap * filt

    if signalTransform != None:
        u *= signalTransform
        cov = filt*numpy.abs(signalTransform)**2
        integral = cov.sum()/mp.Nx/mp.Ny
        u /= integral

    filtered = numpy.real(liteMap.fftTools.ifft2(u))

    filtered = liteMap.liteMapFromDataAndWCS(filtered, mp.wcs)

    if signalTransform == None or forceGaussian == True:
        filtered = filtered.convolveWithGaussian(fwhm, 3)

    return filtered


def filterHorizontalLines( m ):
    stripeProfile = m.data.mean(axis=1)
    a = numpy.ones(m.data.shape[1])
    n = numpy.outer( stripeProfile, a)
    m.data -= n

def filterHorizontalLinesFourier( mp, hw = 20 ):
    ft = fftTools.fftFromLiteMap(mp)
    a = numpy.ones(ft.Ny)
    lx = numpy.outer(a,ft.lx)
    ft.kMap[numpy.where(abs(lx) < hw)] = 0
    mp.data = numpy.real(numpy.fft.ifft2(ft.kMap))

def makeTemplate(m, wl, ell, maxEll, outputFile = None):
    """
    For a given map (m) return a 2D k-space template from a 1D specification wl
    ell = 2pi * i / deltaX
    (m is not overwritten)
    """

    ell = numpy.array(ell)
    wl  = numpy.array(wl)
    
    ft = fftTools.fftFromLiteMap(m)
    #print "max_lx, max_ly", ft.lx.max(), ft.ly.max()
    #print "m_dx, m_dy", m.pixScaleX, m.pixScaleY
    #print "m_nx, m_ny", m.Nx, m.Ny
    l_f = numpy.floor(ft.modLMap)
    l_c = numpy.ceil(ft.modLMap)
    ft.kMap[:] = 0 
    for i in xrange(numpy.shape(ft.kMap)[0]):
        for j in xrange(numpy.shape(ft.kMap)[1]):
            if l_f[i,j] > maxEll or l_c[i,j] > maxEll:
                continue
            w_lo = wl[l_f[i,j]]
            w_hi = wl[l_c[i,j]]
            trueL = ft.modLMap[i,j]
            w = (w_hi-w_lo)*(trueL - l_f[i,j]) + w_lo
            ft.kMap[i,j] = w

    m = m.copy()
    m.data[:] = 0.
    m.data = abs(ft.kMap)
    if outputFile != None:
        m.writeFits(outputFile, overWrite = True)
    return m

def filterWithFunction( mp, fn, maxEll = 39997):
    """
    Filter map mp with function fn
    fn must accept a wavenumber k and return the filter value
    Currently assumes symmetric filter
    """
    ells = numpy.arange(maxEll)
    filt = numpy.zeros(len(ells))
    for i in xrange(len(ells)):
        filt[i] = fn(ells[i])
    print filt
    filt2D = makeTemplate(mp, filt, ells, maxEll)
    ft = numpy.fft.fft2(mp.data)
    ft *= filt2D.data
    mp.data[:] = numpy.real(numpy.fft.ifft2(ft))[:]

def templateFromFourierProfile(mp, model, params, maxEll=39997):
    ell = numpy.arange(maxEll+1)
    wl = model(ell, **params)
    return makeTemplate(mp,wl,ell,maxEll)

def templateFromProfile(mp, model, params, maxEll=39997):
    """
    Filter the map mp with a model profile which accepts
    arguments in degrees in addition to params.
    """
    theta = numpy.linspace(-numpy.pi,numpy.pi,maxEll*2.,endpoint=False)
    m = model(theta*180./numpy.pi, **params)

    x = numpy.fft.fftfreq(mp.Nx)*mp.pixScaleX*mp.Nx
    y = numpy.fft.fftfreq(mp.Ny)*mp.pixScaleY*mp.Ny
    
    dist = ((x[numpy.newaxis,:])**2. + (y[:,numpy.newaxis])**2.)**0.5
    prof = numpy.interp(dist,theta,m)

    out = mp.copy()
    out.data = numpy.real(numpy.fft.fft2(prof))
    
    return out

def estimateBG(dataMap, center_ra, center_dec, radius):
    """ 
    Estimate the background (outside of some disk) as a 2D quadratic function
    @param center_ra deg
    @param center_dec deg
    @param radius deg
    @return background map
    """
    
    d = numpy.copy(dataMap.data)
   
    for j in xrange(dataMap.Nx):
        for i in xrange(dataMap.Ny):
            ra, dec = dataMap.pixToSky(j,i)
            cosdec = numpy.cos(dec*numpy.pi/180.)
            dra =  (ra - center_ra)*cosdec
            ddec = (dec-center_dec)
            dist = (dra**2 + ddec**2)**0.5 
            if dist * 60. < radius:
                d[i,j] = 0
                
    i,j = numpy.where(d != 0.0)
    d2 = d[i,j]
    
    y1, x1 = numpy.shape(d)
    x = numpy.arange(x1) - x1/2
    y = numpy.arange(y1) - y1/2
    xx, yy = numpy.meshgrid(x,y)
    one = numpy.ones(len(i))
    xx2 = xx[i,j]
    yy2 = yy[i,j]

    # calculate coefficients of plane
    M = numpy.array([xx2, yy2, xx2**2, yy2**2, xx2*yy2, one])
    MT = numpy.transpose(M)
    MMT = numpy.dot(M, MT)
    MMT = numpy.linalg.inv(MMT)
    Md = numpy.dot(M, d2)
    K = numpy.dot(MMT, Md)
   
    bg = xx*K[0] + yy*K[1] + +xx*xx*K[2] + yy*yy*K[3] + xx*yy*K[4] + K[5]
    
    bgMap = dataMap.copy()
    bgMap.data[:] = bg[:]

    return bgMap
