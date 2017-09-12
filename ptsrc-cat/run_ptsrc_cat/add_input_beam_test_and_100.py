#This script will be used to add together various fits files for use in source subtration

from astLib import astWCS
from flipper import *
import numpy as np
'''
inFile1 = sys.argv[1]
inFile2 = sys.argv[2]
inFile3 = sys.argv[3]

a = liteMap.liteMapFromFits(inFile1)
b = liteMap.liteMapFromFits(inFile2)
c = liteMap.liteMapFromFits(inFile3)

a.data[:]=a.data[:]+b.data[:]+c.data[:]
a.writeFits("deep56_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set0_2pass_added_I.fits")
'''
inFile1 = "deep56_null_perscanpattern_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set3_2pass_input_used_I.fits"
inFile2 = "deep56_null_perscanpattern_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set3_2pass_100_I.fits"
inFile3 = "deep56_null_perscanpattern_S2_c7v5_AR2_night_nomoon_srcsub_srcpoint_noisecut_fixpickupcut_noturn_debuddyperdet_4way_set3_2pass_beam_test.fits"

a = liteMap.liteMapFromFits(inFile1)
b = liteMap.liteMapFromFits(inFile2)
c = liteMap.liteMapFromFits(inFile3)

a.data[:]=a.data[:]+b.data[:]+c.data[:]
a.writeFits("set3_added_I.fits")