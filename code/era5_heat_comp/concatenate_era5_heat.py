#!/usr/bin/env python

"""
This script concatenates all of the ERA5-HEAT 3-hourly data and puts it in yearly files.

Loading in the 3-hourly data is extremely slow so this is necessary!
"""

import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
from climateforcing.utils import mkdir_p
import numpy as np
from tqdm import tqdm

era5heatdir_out = '/gws/pw/j05/cop26_hackathons/bristol/project10/era5-heat_1deg/'
era5heatdir_in = '/work/scratch-nopw/pmcjs/era5-heat_1deg/'

mkdir_p('/gws/pw/j05/cop26_hackathons/bristol/project10/era5-heat_1deg/')

for year in tqdm(range(1985, 2015)):
    cube_era5 = iris.load(era5heatdir_in + 'ECMWF_utci_%4d*_v1.0_con.nc' % year)
    equalise_attributes(cube_era5)
    unify_time_units(cube_era5)
    for cu in cube_era5:
        cu.coord('time').points = cu.coord('time').points.astype(int)
    cube_era5 = cube_era5.concatenate_cube()
    iris.save(cube_era5, era5heatdir_out + 'ECMWF_utci_%4d_v1.0_con.nc' % year)
