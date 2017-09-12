from flipper import *
import numpy
import cosmology
import scipy.integrate
import scipy.special
from physical_constants import *
from scipy.optimize import fmin, brentq

# A module for creating profiles for, e.g., filtering
def beam_file(x,filename,rescale=1.):
    x = abs(x) # Assumes symmetric
    a,prof = numpy.loadtxt(filename,unpack=True)
    out =  numpy.interp(x,a*rescale,prof,left=0,right=0)
    return out/out.max()

def beta(theta, thetaC = 1./60., beta = 2./3., cutoffStart = 5., cutoffStop = 10.):
    """
    @brief projected isothermal beta model function 
    @param theta angular distance from center in degrees
    @param cutoffCenter smoothly cut off beta at this radius (integral multiple of thetaC )
    @return projected beta model value (center normalized to 1)
    """
    #Build cutoff
    factor = numpy.copy(theta)
    factor[:] = 1.
    cutoffStart *= thetaC
    cutoffStop  *= thetaC
    cutoffWidth = cutoffStop-cutoffStart
    indsCutoff = numpy.where( abs(theta) > cutoffStart )
    factor[indsCutoff] = 0.5*(1+numpy.cos(numpy.pi*(abs(theta[indsCutoff])-cutoffStart)/cutoffWidth))
    inds0 = numpy.where( abs(theta) > cutoffStop)
    factor[inds0] = 0.

    return factor * (1+(theta/thetaC)**2)**(-(3.*beta-1.)/2.)

def gnfw(theta, thetaFWHM = 1./60., alpha = 1.0510, beta = 5.4905, gamma = 0.3081, c500= 1.177, tol = 1e-4, rb=5.):
    def func(x):
        x = numpy.array(x)
        if x.ndim == 0:
            return func([x])[0]
        out = x**-gamma*(1+x**alpha)**((gamma-beta)/alpha)
        out[x==0] = 0.
        return out

    def integrated(b,x_cut=rb*c500):
        if (b > x_cut):
            return b*0.
        x_lim = numpy.sqrt(x_cut**2 - b**2)
        return scipy.integrate.quad(lambda x: func((x**2 + b**2)**.5),
                                0,x_lim,epsrel=tol)[0]
    
    def temp_func(b):
        """Will replace this with a lambda when it stops breaking things"""
        return 0.5 - integrated(b)/integrated(0)

    numpy.seterr(divide='ignore',invalid='ignore') # Get rid of stdout spam
    hwhm = brentq(temp_func, 0, 10)
    fwhm_fac = hwhm*(2./thetaFWHM)
    final = numpy.array([integrated(abs(fwhm_fac*th)) for th in theta])
    numpy.seterr(divide='print',invalid='print') # Restore defaults
    return final/final.max()

def createModelProfile( mp, model, params, ra0 = None, dec0 = None, radius = 0.5 ):
    """
    @brief add a model in a map mp
    @param m map
    @param radius cut off profile at this radius (deg)
    @param model function takes an array of angles (deg) and returns an array of model values
    """

    degToRad = numpy.pi/180.
    if ra0 == None or dec0 == None:
        x0 = mp.Nx/2
        y0 = mp.Ny/2
        ra0,dec0 = mp.pixToSky(x0,y0)
    cosDec0 = numpy.cos(dec0*degToRad)

    ras,decs = liteMap.getCoordinateArrays(mp)
    thetas = (((ras-ra0)*cosDec0)**2 + (decs-dec0)**2)**0.5
    inds = numpy.where(thetas < radius)
    mp.data[inds] += model(thetas[inds], **params)

def generate1DProfileFT( model, params, maxEll = 40000 ):
    """
    Return a profile in fourier space as a function of integer \ell 
    """
    dtheta = 0.5*(2*numpy.pi/maxEll)
    theta = numpy.arange(-numpy.pi, numpy.pi, dtheta)
    m = model(theta*180./numpy.pi, **params)
    ft = abs(numpy.fft.fft(m))
    ell = 2*numpy.pi*numpy.fft.fftfreq(len(ft), dtheta)
    return theta, m, ell, ft


#### For Isothermal Beta-model Case -- Work with Colin Hill ####

f_gas = 0.12 
mu_e  = 1.14 # mean molecular weight per electron of ICM
mu    = 0.59 # mean molecular weight of all ICM species

def R_delta(delta, M, z):
    """
    delta - overdensity wrt background
    M_delta - mass @ this overdensity
    return radius in meters
    """
    return (3*M/(4*numpy.pi)/(delta*cosmology.rho_c(z)))**(1./3.)     

def n_e_0(z,delta,M,r_c,beta):
    """
    z - redshift
    delta - overdensity
    M - Mass within R
    r_c = core radius
    beta - beta model exponent
    """ 
    R_over_r_c = R_delta(delta, M, z)/r_c
    def B(x):
        return 4*numpy.pi*(1+x**2)**(-3*beta/2)*x**2
    _B = scipy.integrate.quad(B, 0, R_over_r_c)
    #print 'B', _B
    return M*f_gas/(mu_e * m_p * r_c**3 * _B[0])

def kT_e(M,z,delta):
    """
    M - Mass
    z - redshift
    delta - overdensity wrt critical density
    """
    return cosmology.G * M * mu * m_p / (3*R_delta(delta, M, z))

def deltaT_beta_physical(angle, z=0.3, delta=200., M=1.e15, vpec=500e3, r_c=.2, beta=1. , nu=148.):
    """
    z - redshift
    delta - overdensity wrt critical
    M - mass at this overdensity in M_solar
    vpec - peculiar velocity in km/s
    r_c - beta model core radius in Mpc
    beta - beta model exponent
    nu - frequency of observations in GHz
    angle - angle off center in arcminutes
    """
    r_c *= Mpc
    M *= m_sun
    nu *= 1e9
    vpec *= 1000.
    #print 'r_c', r_c
    X = h_planck * nu / (k_b * T_cmb)
    _X = X / numpy.tanh(X/2)
    _S = X / numpy.sinh(X/2)
    A = X**4 * numpy.exp(X) / (numpy.exp(X) - 1)**2
    theta = kT_e(M,z,delta)/(m_e*c_light**2)
    #print 'theta', theta
    Y_0 = _X-4.
    #print 'Y_0', Y_0
    Y_1 = -10. + 47./2.*_X - 42./5.*_X**2 + 7./10.*_X**3 + _S**2*(21./5.+7./5.*_X)
    #print 'Y_1', Y_1
    C_0 = 1.
    #print 'C_0', C_0
    C_1 = 10. - 47./5.*_X + 7./5.*_X**2 + 7./10*_S**2
    #print 'C_1', C_1
    #print 'non-rel tsz', theta*Y_0*T_cmb
    #print 'rel tsz', theta*theta*Y_1*T_cmb
    #print 'non-rel ksz', -vpec/c_light*(C_0)*T_cmb
    #print 'rel ksz', -vpec/c_light*(theta*C_1)*T_cmb
    B = (theta*(Y_0+theta*Y_1) - vpec/c_light*(C_0 + theta*C_1))*T_cmb
    tau = (numpy.pi**0.5)*sigma_t*n_e_0(z,delta,M,r_c,beta)*r_c*\
          scipy.special.gamma((3.*beta-1)/2)/scipy.special.gamma(3.*beta/2)
    #print 'ne_0', n_e_0(z,delta,M,r_c,beta)
    #print 'gamma', scipy.special.gamma((3.*beta-1)/2)/scipy.special.gamma(3.*beta/2)
    #print A, B, tau
    da = cosmology.Da(z)
    da *= Mpc

    #print 'angle', angle
    angle_rad = angle*numpy.pi/180./60.
    #print 'rc', r_c/Mpc
    C = (1+da**2*angle_rad**2/r_c**2)**((1-3*beta)/2)
    #print "C", C
    output = A*B*C*tau
    if isinstance(angle,numpy.ndarray):
        inds = numpy.where( abs(angle_rad*da) > R_delta(delta, M, z))
        output[inds] = 0.
    else:
        if abs(angle_rad*da) > R_delta(delta, M, z):
            output = 0.
    return output
