"""
This script is designed to run on the LOTUS cluster on JASMIN.
It will calculate WBGT from the ERA-5 data for a single file. 
It should in theory be applicable this to other models and years, though
caution is needed to make the calendar correct.
For now it saves to /scratch 
Adpated from Chris Smith's UTCI code and scripts

Written by Laila Gohar 31052021
"""
import heatstress_calc
from heatstress_era_calc_tidy import wet_bulb_globe_temperature_era_BoM
import cftime
import datetime
import numpy as np
import iris
import matplotlib.pyplot as pl
import warnings
#import sys

# data output directory: need a GWS!
dataout = "/work/scratch-nopw/lkgohar/"

# data input directory: this will differ for different models
path = "/badc/ecmwf-era5/data/enda/em_sfc/1985/01/01/"


# the final bit of the path will also vary by model and needs to be generalised
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    tas_cube = iris.load_cube(path + 'ecmwf-era5_enda_em_sfc_198501012100*2t.nc')
    tdew_cube = iris.load_cube(path + 'ecmwf-era5_enda_em_sfc_198501012100*2d.nc')


# calculate WBGT
wbgt = wet_bulb_globe_temperature_era_BoM(
    {
        "2t": tas_cube.data,
        "2d": tdew_cube.data,
    }
   
)

# make cubes out of output (only wbgt to save space and time)
wbgt_cube = iris.cube.Cube(
    wbgt,
    units='K',
    var_name='wbgt',
    long_name='Wet Bulb Globe Temperature',
    dim_coords_and_dims=[
        (tas_cube.coord('time'), 0),
        (tas_cube.coord('latitude'), 1),
        (tas_cube.coord('longitude'), 2)
    ]
)

# save the output - should automate by model name etc
iris.save(wbgt_cube, dataout + 'wbgtBoM_3hr_ecmwf-era5_enda_em_sfc_198501012100.nc')
#sys.exit()
