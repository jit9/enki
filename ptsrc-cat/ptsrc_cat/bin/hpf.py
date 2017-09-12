from flipper import liteMap
import numpy
import fff

# mapFile = "map.fits"
# weightFile = "nobs.fits"
# mapFile = "../data/maps2008/diff.fits"
# weightFile = "../data/maps2008/weight.fits"
mapFile = "/scr/queequeg1/shared/maps/temp_map_season_strip_total_cuts3_downweight_0.25_prior_x0_find_10_modes_detau_noise_dedark_noprior_100.fits"
weightFile = "/scr/queequeg1/shared/maps/temp_map_season_strip_total_cuts3_downweight_0.25_prior_x0_find_10_modes_detau_noise_weights.fits"

m = liteMap.liteMapFromFits(mapFile)
w = liteMap.liteMapFromFits(weightFile)
m.data *= numpy.sqrt(w.data/w.data.max())

arcminToRad = 1./60.*numpy.pi/180.
lpf = .75 * arcminToRad
hpf = 2.5 * arcminToRad
# submap = [-30.5, 110.5, -56.7, -48.9]
# submap = [-12, 103, -55.25, -51.1]
submap = None
i = 0
if submap != None:
    m = m.selectSubMap(*submap)
print "pixScaleX", m.pixScaleX
m.data[:] = fff.gaussianFilter(m.data, m.pixScaleX, m.pixScaleY, fwhm=lpf)[:]
m.data[:] = fff.gaussianFilter(m.data, m.pixScaleX, m.pixScaleY, fwhm=hpf, 
        type = 'hpf')[:]
localName = mapFile.split("/")[-1][:-5]
m.writeFits("%s_filt.fits" % localName, overWrite=True)
i+=1

