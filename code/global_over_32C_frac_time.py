from itertools import chain
from glob import glob
import numpy as np
import pandas as pd
import xarray as xr
import os, sys

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

# Mapping libraries:
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

cmap = 'gist_heat_r'
fontsize = 14
cbstr = r'% of record'
projection = 'robinson'
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

    
#project_directory = '/gws/pw/j05/cop26_hackathons/bristol/project10/'
project_directory = os.curdir + '/'

# LAND-SEA MASK:

# I regridded a land-sea mask I have for 20CRv3 to 1x1 and then reset the longitude dimension to match the CMIP6 UTCI dataset format. I use CDO to reset the longitude with:
# cdo sellonlatbox,-180,180,-90,90 landseamask_1x1.nc landseamask.nc

landseamask = xr.load_dataset(project_directory + 'landseamask.nc')
weights = np.cos(np.deg2rad(landseamask.lat))

# Lazy load data

#data_directory = project_directory + 'utci_projections_1deg/BCC-CSM2-MR/historical/r1i1p1f1/'
data_directory = project_directory + ''
filelist = data_directory + 'utci_3hr*.nc'
#paths_to_load = [ glob(f'utci_3hr*.nc') for variable in ['utci'] ]
#paths_to_load = [ glob(filelist) for variable in ['utci'] ]
paths_to_load = glob(filelist)

#dataset = xr.open_mfdataset(paths=chain(*paths_to_load))
#dataset = xr.open_mfdataset(paths_to_load, concat_dim="time", combine="nested", data_vars='minimal', coords='minimal', compat='override', parallel=True)
dataset = xr.open_mfdataset(paths_to_load[0])

# Convert UTCI to degrees Centigrade, set the 32 degree threshold, apply land mask and average over time dimension

utci = dataset.utci[:,:,:]-273.15
threshold = 32.0
utci_over_threshold = xr.where(utci>threshold, utci, np.nan)
utci_over_threshold_frac = (np.isfinite(utci_over_threshold.where(landseamask.LSMASK==1).compute())/len(dataset.time))*100.
utci_over_threshold_frac_mean = utci_over_threshold_frac.mean('time')
utci_over_threshold_frac_weighted = utci_over_threshold_frac.weighted(weights)
utci_over_threshold_frac_weighted_lat = utci_over_threshold_frac_weighted.mean('lon')
utci_over_threshold_frac_weighted_lat_mean = utci_over_threshold_frac_weighted.mean('lon').mean('time')
utci_over_threshold_frac_weighted_mean = utci_over_threshold_frac_weighted.mean(("lon", "lat"))
utci_over_threshold_mean = utci_over_threshold.mean('time').where(landseamask.LSMASK==1).compute()

utci_over_threshold_mean.to_netcdf('global_over_32_mean.nc')
utci_over_threshold_frac.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_mean.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_weighted_lat.to_netcdf('global_over_32_frac.nc')
utci_over_threshold_frac_weighted_mean.to_netcdf('global_over_32_frac.nc')

# PLOTS

fig,ax = plt.subplots(figsize=(15,10))
utci_over_threshold_frac_weighted_mean.plot(label="weighted")
utci_over_threshold_frac.mean(("lon", "lat")).plot(label="unweighted")
plt.ylabel('% of time > 32$^{\circ}$C')
plt.title('BCC-CSM2-MR_historical: 2014-2015 mean UTCI > ' + str(threshold) + r'$^{\circ}C$', fontsize=12)
plt.legend()
plt.savefig('global_over_32C_frac.png')
plt.close('all')

#utci_over_threshold_frac_weighted_mean.where(utci_over_threshold_frac_weighted_mean.time.dt.month==1).plot()
                                
plotfile = 'utci_over_' + str(threshold) + '_' + 'frac' + '.png'
titlestr = 'BCC-CSM2-MR_historical: 2014-2015 mean UTCI > ' + str(threshold) + r'$^{\circ}C$'    

fig, axs = plt.subplots(1,1, figsize=(15,10), subplot_kw=dict(projection=p))
vmin = utci_over_threshold_frac_mean.min(); vmax = utci_over_threshold_frac_mean.max()
g = make_plot(axs, utci_over_threshold_frac_mean, vmin, vmax, cbstr, titlestr, cmap, fontsize)
axs.add_feature(cartopy.feature.OCEAN, zorder=100, alpha=0.2, edgecolor='k')
#cb = fig.colorbar(g, ax=axs.ravel().tolist(), shrink=0.6, extend='both')
cb = fig.colorbar(g, ax=axs, shrink=0.6, extend='both')
cb.set_label(cbstr, rotation=90, labelpad=20, fontsize=fontsize)
cb.ax.tick_params(labelsize=fontsize)
#cb.set_ticks(np.linspace(vmin,vmax,17))
plt.savefig(plotfile, dpi=300)
plt.close('all')

#-----------------------------------------------------------------------------
print('*** END')

print('*** END')


