#uses ndflip to unpack Sigurd's fits files
#saves intensity map in units of Jy/sr and weight map for intensity
import numpy as np
import matplotlib
matplotlib.use('Agg')
from flipper import *
import ndflip

freqGHz=148.65
dir_name="/mnt/act2/mhasse/misc/maps_timestep2/"
map_name="deep5+6+56_tot_sky_map0010.fits"
weight_name="deep5+6+56_tot_sky_div.fits"

maps=ndflip.read_map(dir_name+map_name)
map_I,Q,U=maps
map_I.convertToJyPerSrFromMicroK(freqGHz)

weights =ndflip.read_map(dir_name+weight_name)
wt=weights[0,0]

map_I.writeFits('data_JyPerSr.fits', overWrite=True)
wt.writeFits('weights.fits', overWrite=True)
