import numpy
from physical_constants import *

H_0 = 71.e3 # m/s/Mpc
Omega_m = 0.264
Omega_Lambda =  1- Omega_m
light_speed = 299792458 # m/s

# H^2 = H_0^2 E(z)

def E(z):
    """
    H^2 = H_0^2 E^2(z)
    """
    return numpy.sqrt(Omega_Lambda + Omega_m*(1.+z)**3)

def H(z):
    """
    Hubble Parameter in km/s/Mpc
    """
    return H_0 * E(z)

def Dc(z, dz = 0.0001):
    """
    Comoving Distance (Flat) (Mpc)
    """
    dc = 0.
    zs = numpy.arange(0,z,dz)
    for i in range(len(zs)-1):
        zavg = (zs[i]+zs[i+1])/2
        dc += dz/H(zavg)
    return dc * light_speed

def Da(z, dz = 0.0001):
    """
    Angular Diameter Distance (Da)
    """
    return Dc(z, dz)/(1+z)

def rho_c(z):
    """
    critical density at redshift z
    """
    return 3 * H(z)**2/Mpc**2  / 8 / numpy.pi / G
