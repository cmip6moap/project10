#!/usr/bin/env python

"""
Bias correction for the UTCI dataset

Both the climate model data and the ERA5-HEAT data have been regridded to 1x1 degree and uploaded to JASMIN. Here we use the ERA5-HEAT dataset from 1985 to 2014 and compare this to the derived UTCI from each climate model.

Therefore, instead of bias-correcting temperature or any other variables, we bias correct the derived UTCI.

We therefore assume that ERA5-HEAT is "Truth"! To be fair, I would probably bias correct the individual variables against their ERA5 counterparts. Additionally, for all except temperature this becomes a little tricky and subjective.
"""

import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.pyplot as pl
from climateforcing.utils import mkdir_p
import numpy as np
#import pickle
import scipy.stats as st
from tqdm import tqdm


# ## Obtain historical "training" distributions

era5heatdir = '/gws/pw/j05/cop26_hackathons/bristol/project10/era5-heat_1deg/'
modeldir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/historical/r1i1p1f3/'

## just 30 years for now
## load up the regridding annual chunks and concatenate
cube_era5 = iris.load(era5heatdir + 'ECMWF_utci_*_v1.0_con.nc')
equalise_attributes(cube_era5)
unify_time_units(cube_era5)
for cu in cube_era5:
    cu.coord('time').points = cu.coord('time').points.astype(int)
cube_era5 = cube_era5.concatenate_cube()

## also 30 years of HadGEM3 historical
cube_model = iris.load(modeldir + 'utci_3hr_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_*.nc')
cube_model = cube_model.concatenate_cube()

# generalise this
leeds_model = cube_model[:,143,178]
leeds_era5 = cube_era5[:,143,178]

model_params = {}
model_params['a'] = np.zeros((cube_model.shape[1:3]))
model_params['loc'] = np.zeros((cube_model.shape[1:3]))
model_params['scale'] = np.zeros((cube_model.shape[1:3]))
model_params['lat'] = cube_model.coord('latitude').points
model_params['lon'] = cube_model.coord('longitude').points

era5_params = {}
era5_params['a'] = np.zeros((cube_era5.shape[1:3]))
era5_params['loc'] = np.zeros((cube_era5.shape[1:3]))
era5_params['scale'] = np.zeros((cube_era5.shape[1:3]))
era5_params['lat'] = cube_era5.coord('latitude').points
era5_params['lon'] = cube_era5.coord('longitude').points


model_params['a'][143,178], model_params['loc'][143,178], model_params['scale'][143,178] = st.skewnorm.fit(leeds_model.data)
era5_params['a'][143,178], era5_params['loc'][143,178], era5_params['scale'][143,178] = st.skewnorm.fit(leeds_era5.data)

# ## How to bias correct
# 
# $\hat{x}_{m,p}(t) = F^{-1}_{o,h} ( F_{m,h} (x_{m,p}(t)) )$
# 
# - $x_{m,p}$ is the future predicted variable, i.e. the SSP value from the climate model
# - $F_{m,h}$ is the CDF of the historical period in the climate model
# - $F_{o,h}$ is the CDF of the historical period in the observations (or in this case, ERA5)


# F_{m,h}
# In: st.skewnorm.cdf(290, model_params['a'][143,178], model_params['loc'][143,178], model_params['scale'][143,178])
# Out: 0.4921534798137802   # percentile of 290 K in HadGEM3 climate

# F^{-1}_{o,h}
# In: st.skewnorm.ppf(0.4921534798137802, era5_params['a'][143,178], era5_params['loc'][143,178], era5_params['scale'][143,178])
# Out: 290.57999427509816  # UTCI in ERA5 corresponding to this percentile.

# transfer function
def bias_correct(x, model_params, obs_params, ilat, ilon):
    cdf = st.skewnorm.cdf(x, model_params['a'][ilat, ilon], model_params['loc'][ilat, ilon], model_params['scale'][ilat, ilon])
    x_hat = st.skewnorm.ppf(cdf, obs_params['a'][ilat, ilon], obs_params['loc'][ilat, ilon], obs_params['scale'][ilat, ilon])
    return x_hat


# ## Bias correct future simulations
# 
# For now, just use 2100

modelfuturedir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/ssp585/r1i1p1f3/'
cube_model_future = iris.load(modelfuturedir + 'utci_3hr_HadGEM3-GC31-LL_ssp585_r1i1p1f3_gn_210001010300-210101010000.nc')
cube_model_future = cube_model_future.concatenate_cube()


leeds_model_future = cube_model_future[:,143,178]


model_future_params = {}
model_future_params['a'] = np.zeros((cube_model_future.shape[1:3]))
model_future_params['loc'] = np.zeros((cube_model_future.shape[1:3]))
model_future_params['scale'] = np.zeros((cube_model_future.shape[1:3]))
model_future_params['lat'] = cube_model_future.coord('latitude').points
model_future_params['lon'] = cube_model_future.coord('longitude').points

model_future_params['a'][143,178], model_future_params['loc'][143,178], model_future_params['scale'][143,178] = st.skewnorm.fit(leeds_model_future.data)

#pl.hist(leeds_model.data, density=True, label='HadGEM3-GC31-LL 1985', alpha=0.3, bins=50)
#pl.hist(leeds_era5.data, density=True, label='ERA5-HEAT', alpha=0.3, bins=50)
#pl.hist(leeds_model_future.data, density=True, label='HadGEM3-GC31-LL 2100', alpha=0.3, bins=50)
#pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), model_params['a'][143,178], model_params['loc'][143,178], model_params['scale'][143,178]), color='tab:blue')
#pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), era5_params['a'][143,178], era5_params['loc'][143,178], era5_params['scale'][143,178]), color='tab:orange')
#pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), model_future_params['a'][143,178], model_future_params['loc'][143,178], model_future_params['scale'][143,178]), color='tab:green')
#pl.legend()
#pl.title('Leeds grid cell')
#pl.show()


# bias correct the Leeds 2100 projections
leeds_model_future_biascorrected = bias_correct(leeds_model_future.data, model_params, era5_params, 143, 178)

pl.hist(leeds_model.data, density=True, label='HadGEM3-GC31-LL 1985', alpha=0.3, bins=50)
pl.hist(leeds_era5.data, density=True, label='ERA5-HEAT', alpha=0.3, bins=50)
pl.hist(leeds_model_future.data, density=True, label='HadGEM3-GC31-LL 2100', alpha=0.3, bins=50)
pl.hist(leeds_model_future_biascorrected, density=True, label='Bias-corrected 2100', alpha=0.3, bins=50)

pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), model_params['a'][143,178], model_params['loc'][143,178], model_params['scale'][143,178]), color='tab:blue')
pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), era5_params['a'][143,178], era5_params['loc'][143,178], era5_params['scale'][143,178]), color='tab:orange')
pl.plot(np.arange(240, 320), st.skewnorm.pdf(np.arange(240, 320), model_future_params['a'][143,178], model_future_params['loc'][143,178], model_future_params['scale'][143,178]), color='tab:green')
pl.legend()
pl.title('Leeds grid cell')
pl.show()
