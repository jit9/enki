from flipper import *
import scipy.ndimage

def optimalFilter(mp, signalTransform = None):
    """
    @param mp liteMap to filter
    @param noiseSpectrum an array containing the noise power spectrum of mp
    @param signalTransform an array containing the fourier transform of the signal for matched filtering
    @param fwhm if signalTransform is not specified, then filter with a gaussian of this width
    @param extraFilt another filter array to apply (form = liteMap)
    """

    """ get power spectra of cmb and noise maps """
    pNoise = liteMap.fftTools.powerFromLiteMap(mp)

    """ calculate filter """
    filt = 1 / pNoise.powerMap

#if we are not using a nice model, the filter need conditioning...
#eliminate outliers from filter
    med = numpy.median(filt)
    filt[numpy.where(filt> 10*med)] = med
    # convolve filter with gaussian to smooth it
    kernel_width = (5,5)  # width of gaussian kernal in sigma
    filt = scipy.ndimage.gaussian_filter(filt, kernel_width)

    ft = liteMap.fftFromLiteMap(mp)
    u = ft.kMap * filt

    if signalTransform != None:
        u *= signalTransform
        cov = filt*numpy.abs(signalTransform)**2
        integral = cov.sum()/mp.Nx/mp.Ny
        u /= integral

    filtered = numpy.real(liteMap.fftTools.ifft2(u))

    filtered = liteMap.liteMapFromDataAndWCS(filtered, mp.wcs)

    return filtered

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
