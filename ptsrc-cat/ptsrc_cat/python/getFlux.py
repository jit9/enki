from numpy import *
from scipy.linalg import *
from cluster import *
from flipper import *

def getFlux(dataMap, ra, dec, radius, objectType='source'):
    
    # cut small submap
    #ra = cl.ra
    #dec = cl.dec
    centX, centY = dataMap.skyToPix(ra, dec)
    x0 = centX - 20
    x1 = centX + 20
    y0 = centY - 20
    y1 = centY + 20
    ra1, dec1 = dataMap.pixToSky(x0, y0)
    ra2, dec2 = dataMap.pixToSky(x1, y1)

    subLarge = dataMap.selectSubMap(ra1, ra2, dec1, dec2)

    if objectType == 'source':
        subLarge.convertToJyPerSrFromMicroK(148)

    d = copy(subLarge.data)
    d[10:32, 10:32] = 0
    i, j = where(d != 0.0)
    d2 = d[i,j]
    
    x1, y1 = shape(d)
    x = arange(0., x1) - x1/2
    y = arange(0., y1) - y1/2
    xx, yy = meshgrid(x,y)
    one = ones(len(i))
    xx2 = xx[i,j]
    yy2 = yy[i,j]

    # calculate coefficients of plane
    M = array([xx2, yy2, xx2**2, yy2**2, xx2*yy2, one])
    MT = transpose(M)
    MMT = dot(M, MT)
    MMT = inv(MMT)
    Md = dot(M, d2)
    K = dot(MMT, Md)
   
    plane = xx*K[0] + yy*K[1] + +xx*xx*K[2] + yy*yy*K[3] + xx*yy*K[4] + K[5]

    # subtract off plane
    new = subLarge.data - plane

    mapNew = liteMap.liteMapFromDataAndWCS(new, subLarge.wcs)

    
    # sum up pixels within given radius
    ra1 = ra + radius/60.
    dec1 = dec + radius/60.
    x1, y1 = mapNew.skyToPix(ra1, dec1)
    centX, centY = mapNew.skyToPix(ra, dec)

    xDist = abs(x1 - centX)
    yDist = abs(y1 - centY)
    x_, y_ = shape(d)
    x = arange(0., x_)
    y = arange(0., y_) 
    xx, yy = meshgrid(x, y)
    
    Dist = sqrt((xx-centX)**2 + (yy-centY)**2)
    i, j = where(Dist < yDist)
    
    
    # integrate pixels to get flux
    C = sum(new[i, j])    
    A = C*mapNew.pixScaleX*mapNew.pixScaleY*(180/pi*60)**2
    B = C*mapNew.pixScaleX*mapNew.pixScaleY

    # get frequency dependence
    T_cmb = 2.725
    k = 1.3806503e-23
    h = 6.626068e-34

    nu = 148 #in GHz
    x = (h*nu)/(k*T_cmb)*1e9
    f = x*(numpy.exp(x) + 1)/(numpy.exp(x) - 1) - 4

    # y1 in arcmin**2 and y2 in steradians
    y1 = A/((2.725*10**6)*f)
    #y2 = B/((2.725*10**6)*f)

    if objectType == 'source':
        return B
    elif objectType == 'cluster':
        return y1
    
