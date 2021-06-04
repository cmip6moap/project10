import cartopy.crs as ccrs
import matplotlib.pyplot as plt

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib import colors


def plot_diff_maps_monthly(data, title, vmin=-10, vmax=10, cmap_diverging=True, outputname=None):
    '''
    Plot pcolormesh maps for each month as 12 panels. Expects input xarray.DataArray with "month"
    dimension as the time co-ordinate e.g.
    
        data (month, lat, lon)
        
    Args:
        data (xarray.DataArray) :
            DataArray or single entry selected from a Dataset. This should contain spatial data for
            each month.
        title (str) :
            Title to include on the plot
        vmin, vmax (float/int, optional) :
            Minimum and maximum values to display. Colourbar is set to extend on both ends, so any values
            outside of the threshold will be shown as the max or min values in the colormap.
            Note if cmap_diverging is set to True, vmin must be < 0.0 and vmax must be > 0.0
            Default = -10, 10
        cmap_diverging (bool, optional) :
            Use Diverging colourmap - i.e. a colourmap with white in the middle. To allow both positive and
            negative values set this to True. Otherwise, set this to False.
            The diverging colourmap used is "coolwarm"
            The non-diverging colourmap used is "inferno_r"
            Default = True
        outputname (str/None, optional) :
            To write the output to file set the outputname. This should include the extension e.g. .png
            If this is None, no output will be written and will only be interactively displayed.
        
    Returns:
        None
        
        Plots to screen
        If outputname is specified, also saves to file.
    '''
    
    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(14,12), subplot_kw={"projection":ccrs.PlateCarree()})
    axes = axes.flatten()

    months = data.month
    # month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    if cmap_diverging:
        offset = colors.TwoSlopeNorm(vmin=vmin, vcenter=0., vmax=vmax)
    
    # Used suggestions from http://xarray.pydata.org/en/stable/examples/monthly-means.html
    for i, month in enumerate(months):
        
        if cmap_diverging:
            data.sel(month=month).plot.pcolormesh(ax=axes[i], cmap="coolwarm", 
                                                  norm=offset, extend="both",
                                                  add_colorbar=True)
        else:
            data.sel(month=month).plot.pcolormesh(ax=axes[i], cmap="inferno_r", 
                                                  vmin=vmin, vmax=vmax, extend="both",
                                                  add_colorbar=True)

        axes[i].coastlines() # Show edges of land
        axes[i].add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k', facecolor='0.9') # Mask the ocean
        axes[i].set_extent([-179.5, 180, -68, 90]) # Cut off Antartica
        #axes[i].set_ylabel(month_names[i]) # Wasn't displaying (or was displaying underneath something else)

    fig.tight_layout()
    fig.suptitle(title, fontsize=16, y=1.02)
    
    fig.show()
    
    if outputname:
        fig.savefig(outputname)