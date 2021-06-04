#!/usr/bin/env python

"""
This script regrids the climate model output data to 1x1 degree.
"""

import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import iris.analysis.cartography
import iris.coord_categorisation
import numpy as np
import glob
from climateforcing.utils import mkdir_p
import sys

# Make a dummy grid
latitude = iris.coords.DimCoord(
    np.arange(-89.5,90,1),
    standard_name='latitude',
    units='degrees',
    long_name='Latitude',
    var_name='lat',
    coord_system=None
)
longitude = iris.coords.DimCoord(
    np.arange(-179.5,180,1),
    standard_name='longitude',
    long_name='Longitude',
    var_name='lon',
    units='degrees',
    circular=True,
    coord_system=None
)

ny = len(latitude.points)
nx = len(longitude.points)

dummy_data = np.zeros((ny, nx))
dummy_cube = iris.cube.Cube(dummy_data, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
dummy_cube.coord('longitude').guess_bounds()
dummy_cube.coord('latitude').guess_bounds()

# do the regrid and save output
origdir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections/'
regriddir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/'

models = ['HadGEM3-GC31-LL', 'BCC-CSM2-MR']
runs = {
    'HadGEM3-GC31-LL': ['r1i1p1f3'],
    'BCC-CSM2-MR': ['r1i1p1f1'],
}
scenarios = {}
scenarios['HadGEM3-GC31-LL'] = {}
scenarios['BCC-CSM2-MR'] = {}
scenarios['HadGEM3-GC31-LL']['r1i1p1f3'] = ['historical', 'ssp126', 'ssp245', 'ssp585']
scenarios['BCC-CSM2-MR']['r1i1p1f1'] = ['historical', 'ssp126', 'ssp245', 'ssp585']

for model in models:
    for run in runs[model]:
        first_file=True
        for scenario in scenarios[model][run]:
            filelist = glob.glob(origdir + '%s/%s/%s/*.nc' % (model, scenario, run))
            for file in filelist:
                cube_orig = iris.load_cube(file)
                filename = file.split('/')[-1]
                if not cube_orig.coord('longitude').has_bounds():
                    cube_orig.coord('longitude').guess_bounds()
                if not cube_orig.coord('latitude').has_bounds():
                    cube_orig.coord('latitude').guess_bounds()

                # regrid the model to the new grid
                if first_file:
                    regridder = iris.analysis.AreaWeighted().regridder(cube_orig, dummy_cube)
                    first_file = False
                cube_regrid = regridder(cube_orig)
                # save the output
                outdir = regriddir + '%s/%s/%s/' % (model, scenario, run)
                mkdir_p(outdir)
                iris.save(cube_regrid, outdir + filename)
                sys.stdout.write(filename + ' success\n') 
