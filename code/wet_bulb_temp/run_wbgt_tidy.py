#!/usr/bin/env python
"""
This script is designed to run on the LOTUS cluster on JASMIN.

It will calculate WBGT from a specified model and under a specified scenario passed
on the command line.

For now it saves to /scratch (we need a group workspace for this work)

Usage:
    python run_wbgt.py index

index is an integer that provides the model, scenario and run to use and is a lookup to
a csv file.

model: precise name of the model to run (e.g IITM-ESM)
scenario: precise CMIP6-style name of the scenario to run (e.g. ssp585)
run: the realisation of the scenario (i.e. r1i1p1f3)

This is adapted from Chris Smith's excellant UTCI code and scripts.

Written by Laila Gohar 25052021 laila.gohar@metoffice.gov.uk
"""

#from climateforcing.utci import universal_thermal_climate_index, mean_radiant_temperature

from climateforcing.solar import cos_mean_solar_zenith_angle, modified_julian_date
from climateforcing.utils import mkdir_p
import cftime
import datetime
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import matplotlib.pyplot as pl
import pandas as pd
import warnings
import sys
import glob
import heatstress_calc
from heatstress_calc_210521 import wet_bulb_globe_temperature_BoM

index = sys.argv[1]

df = pd.read_csv('3hr_models.csv', index_col=0)
model, scenario, run, startyear, endyear = df.loc[int(index)]


print(model, scenario, run)
print()

# data output directory: need a GWS!
dataout = "/work/scratch-nopw/lkgohar/%s/%s/%s/" % (model, scenario, run)
mkdir_p(dataout)

# data input directory: this will differ for different models. With iris we can use
# glob
path = "/badc/cmip6/data/CMIP6/*/*/%s/%s/%s/3hr/" % (model, scenario, run)

# the final bit of the path will also vary by model grid
vars = ['tas', 'huss']
cubes = {}
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for var in vars:
        cubes[var] = iris.load(path + '%s/*/latest/%s_3hr_%s_%s_%s_*_*.nc' % (var, var, model, scenario, run))
        equalise_attributes(cubes[var])
        unify_time_units(cubes[var])
        cubes[var] = cubes[var].concatenate_cube()
fp = glob.glob(path + '%s/*/latest/%s_3hr_%s_%s_%s_*_*.nc' % (var, var, model, scenario, run))
grid = fp[0].split('/')[12]

# check that all of the cubes are the same number of time points.
# If they are not, there are some missing variable slices and
# the outputs will not make sense.
# it's sufficient to check everything relative to tas
tas_shape = cubes['tas'].shape[0]
for var in vars:
    if cubes[var].shape[0] != tas_shape:
        raise ValueError(var + ' is a different number of time points to tas.')

# get lat and lon
lat = cubes['tas'].coord('latitude').points
lon = cubes['tas'].coord('longitude').points
lonmesh, latmesh = np.meshgrid(lon, lat)

nlat = len(lat)
nlon = len(lon)

# sort out times
all_time_coord = cubes['tas'].coord('time')
all_time_points = all_time_coord.units.num2date(all_time_coord.points)
n_all_time = len(all_time_points)
calendar = all_time_coord.units.calendar
first_time = all_time_points[0]
last_time = all_time_points[-1]
#sys.exit()
if scenario == 'historical':
   if last_time.year < 2013: sys.exit()
else:
   if last_time.year < 2098: sys.exit()
# start the year chunking loop
for year in range(int(startyear), int(endyear)):
    # historical: we want to start in 1985
#    if year < 1985:
#        continue
    # tas timesteps are usually at the end of the period, i.e. 03:00 for 00:00 to 03:00.
    # the first timestep of each year is thus 01 January at 03:00 UTC.
    # it is not clear whether this is an 03:00 instantaneous value or 00:00 to 03:00 mean.
    # hopefully, if it is instantaneous, 3hr data means the biases won't be huge. It also
    # means that the radiation and temperature time steps are in sync.

    # sometimes however this isn't the case. The BCC model has tas timesteps running from
    # 00:00 which are 22:30 to 01:30 means. In BCC the first time step of the year is 
    # 01 January at 00:00 UTC. The MRI model also starts at 00:00 rather than 03:00, and
    # the metadata confirms that these are point values and not time means in MRI. It
    # also means that in these two models the temperature and radiation are 90 mins out
    # of sync - hopefully this is not critical, but if we wanted to adjust for this we
    # could take the mean of two neighbouring time points.
    # there's probably a nice iris-friendly way to do this, but for now, just treat BCC
    # and MRI as exceptions to the usual rule, and be sure to check the filenames and
    # metadata for any new models.
    if model in ['BCC-CSM2-MR', 'MRI-ESM2-0']:
        first_hour_tas=0
    else:
        first_hour_tas=3

    i_start = int(8 * (cftime.date2num(cftime.datetime(year, 1, 1, first_hour_tas, 0, 0, calendar=calendar), 'days since 1850-01-01 00:00', calendar=calendar) - cftime.date2num(first_time, 'days since 1850-01-01 00:00', calendar=calendar)))
    print(model, scenario, run, year, i_start)
    days_in_year = cftime.date2num(cftime.datetime(year+1, 1, 1, 0, 0, 0, calendar=calendar), 'days since 1850-01-01 00:00', calendar=calendar) - cftime.date2num(cftime.datetime(year, 1, 1, 0, 0, 0, calendar=calendar), 'days since 1850-01-01 00:00', calendar=calendar)
    timepoints_in_year = 8 * days_in_year  # 8 x 3-hr timepoints per day
    i_end = i_start + timepoints_in_year
    mjd = np.zeros(timepoints_in_year)
    # radiation timepoints all seem to be 01:30, 04:30, ...
    first_day_of_year = cftime.datetime(year, 1, 1, 1, 30, calendar=calendar)
    for imjd in range(timepoints_in_year):
        mjd[imjd] = modified_julian_date(first_day_of_year + datetime.timedelta(hours=3*imjd))


    # calculate WBGT
    wbgt = wet_bulb_globe_temperature_BoM(
        {
            "tas": cubes['tas'][i_start:i_end,...].data,
            "huss": cubes['huss'][i_start:i_end,...].data,
            "ps": cubes['ps'][i_start:i_end,...].data
        }
    )

    # make cubes out of output (only WBGT to save space and time)
    wbgt_cube = iris.cube.Cube(
        wbgt,
        units='K',
        var_name='wbgt',
        long_name='Wet Bulb Globe Temperature',
        dim_coords_and_dims=[
            (cubes['tas'][i_start:i_end,...].coord('time'), 0),
            (cubes['tas'][i_start:i_end,...].coord('latitude'), 1),
            (cubes['tas'][i_start:i_end,...].coord('longitude'), 2)
        ]
    )
    # save the output
    iris.save(wbgt_cube, dataout + 'wbgt_3hr_%s_%s_%s_%s_%4d01010300-%d01010000.nc' % (model, scenario, run, grid, year, year+1))
