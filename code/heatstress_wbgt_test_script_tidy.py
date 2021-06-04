"""
This script is designed to run on the LOTUS cluster on JASMIN.
It will calculate WBGT from the HadGEM3-GC31-LL model for SSP5-8.5 for
the year 2100. It should in theory be applicable this to other models and years, though
caution is needed to make the calendar correct.
For now it saves to /scratch. Adapted from Chris Smith's script for WBGT test calculation

"""
import heatstress_calc_tidy
from heatstress_calc_tidy import wet_bulb_globe_temperature_BoM
import cftime
import datetime
import numpy as np
import iris
import matplotlib.pyplot as pl
import warnings

# data output directory: need a GWS!
dataout = "/work/scratch-nopw/lkgohar/"

# data input directory: this will differ for different models
path = "/badc/cmip6/data/CMIP6/ScenarioMIP/MOHC/HadGEM3-GC31-LL/ssp585/r1i1p1f3/3hr/"


# Read in input data 
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    tas_cube = iris.load_cube(path + 'tas/gn/latest/tas_3hr_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_210001010300-210101010000.nc')
    huss_cube = iris.load_cube(path + 'huss/gn/latest/huss_3hr_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_210001010300-210101010000.nc')
    ps_cube = iris.load_cube(path + 'ps/gn/latest/ps_3hr_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_210001010300-210101010000.nc')


# calculate WBGT
wbgt = wet_bulb_globe_temperature_BoM(
    {
        "tas": tas_cube.data,
        "huss": huss_cube.data,
	"ps": ps_cube.data
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
iris.save(wbgt_cube, dataout + 'wbgtBoM_3hr_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_210001010300-210101010000.nc')
