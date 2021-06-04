#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: utci_exceedence_anomaly.py
#------------------------------------------------------------------------------
# Version 0.1
# 4 June, 2021
# Michael Taylor
# https://patternizer.github.io
# michael DOT a DOT taylor AT uea DOT ac DOT uk 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------

from itertools import chain
from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
import os, sys
from pathlib import Path

# Plotting libraries:
import matplotlib
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import colors as mplc
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.ticker as mticker

# %matplotlib inline

# Mapping libraries:
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Geo libraries:

import geopandas as gp
#import pooch
#import pyproj
#from shapely.geometry import Polygon
#conda install -c conda-forge regionmask
#pip install git+https://github.com/mathause/regionmask
import regionmask

from scipy import stats

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

#----------------------------------------------------------------------------
# SETTINGS
#----------------------------------------------------------------------------

cmap = 'coolwarm'
fontsize = 14
cbstr = r'UTCI [$^{\circ}$C]'
flag_ar6 = False # False -->  use AR5
projection = 'equalearth'
if projection == 'platecarree': p = ccrs.PlateCarree(central_longitude=0); threshold = 0
if projection == 'mollweide': p = ccrs.Mollweide(central_longitude=0); threshold = 1e6
if projection == 'robinson': p = ccrs.Robinson(central_longitude=0); threshold = 0
if projection == 'equalearth': p = ccrs.EqualEarth(central_longitude=0); threshold = 0
if projection == 'geostationary': p = ccrs.Geostationary(central_longitude=0); threshold = 0
if projection == 'goodehomolosine': p = ccrs.InterruptedGoodeHomolosine(central_longitude=0); threshold = 0
if projection == 'europp': p = ccrs.EuroPP(); threshold = 0
if projection == 'northpolarstereo': p = ccrs.NorthPolarStereo(); threshold = 0
if projection == 'southpolarstereo': p = ccrs.SouthPolarStereo(); threshold = 0
if projection == 'lambertconformal': p = ccrs.LambertConformal(central_longitude=0); threshold = 0

#------------------------------------------------------------------------------
# METHODS
#------------------------------------------------------------------------------
    
# Xarray plotting function with 5x5 grid overlay

def make_plot(axi,v,vmin,vmax,cbstr,titlestr,cmap,fontsize):
    
    g = v.plot(ax=axi, transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=cmap, cbar_kwargs={'orientation':'vertical','extend':'both','shrink':0.7, 'pad':0.1})         
    cb = g.colorbar; cb.ax.tick_params(labelsize=fontsize); cb.set_label(label=cbstr, size=fontsize); cb.remove()
    axi.set_global()        
    axi.coastlines(color='grey')
    axi.set_title(titlestr, fontsize=fontsize)    
    gl = axi.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, linewidth=1, color='black', alpha=0.2, linestyle='-')
    gl.top_labels = False; gl.bottom_labels = False; gl.left_ylabels = False; gl.right_ylabels = False
    gl.xlines = True; gl.ylines = True
    gl.xlocator = mticker.FixedLocator(np.linspace(-180,180,73)) # every 5 degrees
    gl.ylocator = mticker.FixedLocator(np.linspace(-90,90,37))   # every 5 degrees
    gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER

    return g

#------------------------------------------------------------------------------
# SET: paths
#------------------------------------------------------------------------------

#project_directory = '/gws/pw/j05/cop26_hackathons/bristol/project10/'
project_directory = os.curdir + '/'

#data_directory = project_directory + 'utci_projections_1deg/BCC-CSM2-MR/historical/r1i1p1f1/'
data_directory = project_directory + 'MONTHLY/utci_projections_1deg_monthly/'

#------------------------------------------------------------------------------
# LOAD: landseamask & latitudinal weights
#------------------------------------------------------------------------------

landseamask = xr.load_dataset(project_directory + 'landseamask.nc')
weights = np.cos(np.deg2rad(landseamask.lat))

#------------------------------------------------------------------------------
# CALCULATIONS
#------------------------------------------------------------------------------

# Convert UTCI to degrees Centigrade, set the 32 degree threshold, apply land mask and average over time dimension

threshold = 32.0

baseline = xr.open_dataset(data_directory + '/HadGEM3-GC31-LL/historical/r1i1p1f3/monthly_avg.nc')
baseline_land = baseline.utci.where(landseamask.LSMASK==1)-273.15
baseline_land_over_threshold = xr.where(baseline_land>threshold, baseline_land, np.nan)
baseline_normals = baseline_land_over_threshold.sel(time=slice("1986", "2016")).groupby('time.month').mean(['time'])
                                          
ssp126 = xr.open_dataset(data_directory + '/HadGEM3-GC31-LL/ssp126/r1i1p1f3/monthly_avg.nc')
ssp245 = xr.open_dataset(data_directory + '/HadGEM3-GC31-LL/ssp245/r1i1p1f3/monthly_avg.nc')
ssp585 = xr.open_dataset(data_directory + '/HadGEM3-GC31-LL/ssp585/r1i1p1f3/monthly_avg.nc')

ssp126_land = ssp126.utci.where(landseamask.LSMASK==1)-273.15
ssp245_land = ssp245.utci.where(landseamask.LSMASK==1)-273.15
ssp585_land = ssp585.utci.where(landseamask.LSMASK==1)-273.15

ssp126_land_over_threshold = xr.where(ssp126_land>threshold, ssp126_land, np.nan)
ssp245_land_over_threshold = xr.where(ssp245_land>threshold, ssp245_land, np.nan)
ssp585_land_over_threshold = xr.where(ssp585_land>threshold, ssp585_land, np.nan)

ssp126_utci = ssp126_land_over_threshold.groupby('time.month') - baseline_normals
ssp245_utci = ssp245_land_over_threshold.groupby('time.month') - baseline_normals
ssp585_utci = ssp585_land_over_threshold.groupby('time.month') - baseline_normals

ssp126_utci_over_threshold_mean = ssp126_utci.mean('time')
ssp245_utci_over_threshold_mean = ssp245_utci.mean('time')
ssp585_utci_over_threshold_mean = ssp585_utci.mean('time')

#------------------------------------------------------------------------------
# PLOTS
#------------------------------------------------------------------------------

#plotfile = 'bcc' + '_' + 'ssp126' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#plotfile = 'bcc' + '_' + 'ssp245' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#plotfile = 'bcc' + '_' + 'ssp585' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#plotfile = 'cmcc' + '_' + 'ssp126' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#plotfile = 'cmcc' + '_' + 'ssp245' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#plotfile = 'cmcc' + '_' + 'ssp585' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
#titlestr = 'BCC-CSM2-MR: SSP126 (2015-2100) UTCI>32$^{\circ}$C anomaly'    
#titlestr = 'BCC-CSM2-MR: SSP245 (2015-2100) UTCI>32$^{\circ}$C anomaly'    
#titlestr = 'BCC-CSM2-MR: SSP585 (2015-2100) UTCI>32$^{\circ}$C anomaly'    
#titlestr = 'CMCC-ESM2: SSP126 (2015-2100) UTCI>32$^{\circ}$C anomaly'    
#titlestr = 'CMCC-ESM2: SSP245 (2015-2100) UTCI>32$^{\circ}$C anomaly'    
#titlestr = 'CMCC-ESM2: SSP585 (2015-2100) UTCI>32$^{\circ}$C anomaly'    

plotfile = 'hadgem3' + '_' + 'ssp126' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
titlestr = 'HadGEM3-GC31-LL: SSP126 (2015-2100) UTCI>32$^{\circ}$C anomaly'    

fig, axs = plt.subplots(1,1, figsize=(15,10), subplot_kw=dict(projection=p))
#vmin = np.nanmin(ssp126_utci_over_threshold_mean); vmax = np.nanmax(ssp126_utci_over_threshold_mean)
vmin=-5; vmax=5
g = make_plot(axs, ssp126_utci_over_threshold_mean, vmin, vmax, cbstr, titlestr, cmap, fontsize)
axs.add_feature(cartopy.feature.OCEAN, zorder=100, alpha=0.2, edgecolor='k')
#cb = fig.colorbar(g, ax=axs.ravel().tolist(), shrink=0.6, extend='both')
cb = fig.colorbar(g, ax=axs, shrink=0.6, extend='both')
cb.set_label(cbstr, rotation=90, labelpad=20, fontsize=fontsize)
cb.ax.tick_params(labelsize=fontsize)
#cb.set_ticks(np.linspace(vmin,vmax,17))
plt.savefig(plotfile, dpi=300)
plt.close('all')

plotfile = 'hadgem3' + '_' + 'ssp245' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
titlestr = 'HadGEM3-GC31-LL: SSP245 (2015-2100) UTCI>32$^{\circ}$C anomaly'    

fig, axs = plt.subplots(1,1, figsize=(15,10), subplot_kw=dict(projection=p))
#vmin = np.nanmin(ssp245_utci_over_threshold_mean); vmax = np.nanmax(ssp245_utci_over_threshold_mean)
vmin=-5; vmax=5
g = make_plot(axs, ssp245_utci_over_threshold_mean, vmin, vmax, cbstr, titlestr, cmap, fontsize)
axs.add_feature(cartopy.feature.OCEAN, zorder=100, alpha=0.2, edgecolor='k')
#cb = fig.colorbar(g, ax=axs.ravel().tolist(), shrink=0.6, extend='both')
cb = fig.colorbar(g, ax=axs, shrink=0.6, extend='both')
cb.set_label(cbstr, rotation=90, labelpad=20, fontsize=fontsize)
cb.ax.tick_params(labelsize=fontsize)
#cb.set_ticks(np.linspace(vmin,vmax,17))
plt.savefig(plotfile, dpi=300)
plt.close('all')

plotfile = 'hadgem3' + '_' + 'ssp585' + '_' + 'utci_over_' + str(threshold) + 'anomaly' '.png'
titlestr = 'HadGEM3-GC31-LL: SSP585 (2015-2100) UTCI>32$^{\circ}$C anomaly'    

fig, axs = plt.subplots(1,1, figsize=(15,10), subplot_kw=dict(projection=p))
#vmin = np.nanmin(ssp585_utci_over_threshold_mean); vmax = np.nanmax(ssp585_utci_over_threshold_mean)
vmin=-5; vmax=5
g = make_plot(axs, ssp585_utci_over_threshold_mean, vmin, vmax, cbstr, titlestr, cmap, fontsize)
axs.add_feature(cartopy.feature.OCEAN, zorder=100, alpha=0.2, edgecolor='k')
#cb = fig.colorbar(g, ax=axs.ravel().tolist(), shrink=0.6, extend='both')
cb = fig.colorbar(g, ax=axs, shrink=0.6, extend='both')
cb.set_label(cbstr, rotation=90, labelpad=20, fontsize=fontsize)
cb.ax.tick_params(labelsize=fontsize)
#cb.set_ticks(np.linspace(vmin,vmax,17))
plt.savefig(plotfile, dpi=300)
plt.close('all')

#------------------------------------------------------------------------------
# WORK IN PROGRESS
#------------------------------------------------------------------------------

# Extract fraction of time UTCI is above threshold, weight and slice by latitude

ssp126_utci_over_threshold_frac = (np.isfinite(ssp126_utci_over_threshold_mean)/len(ssp126_utci.time))*100.
utci_over_threshold_frac_mean = utci_over_threshold_frac.mean('time') # time-averaged map
utci_over_threshold_frac_weighted = utci_over_threshold_frac.weighted(weights)
utci_over_threshold_frac_weighted_lat = utci_over_threshold_frac_weighted.mean('lon')
utci_over_threshold_frac_weighted_lat_mean = utci_over_threshold_frac_weighted.mean('lon').mean('time')
utci_over_threshold_frac_weighted_mean = utci_over_threshold_frac_weighted.mean(("lon", "lat"))

# SAVE: extracts to netCDF

utci_over_threshold_mean.to_netcdf('global_over_32_mean.nc')
utci_over_threshold_frac.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_mean.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_weighted_lat.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_weighted_mean.to_netcdf('global_over_32_frac.nc')

#-----------------------------------------------------------------------------
print('*** END')


