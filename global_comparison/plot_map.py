import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def plot_diff_maps_monthly(data, title, vmin=-10, vmax=10, ):

    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(14,12), subplot_kw={"projection":ccrs.PlateCarree()})
    axes = axes.flatten()

    months = data.month
    #month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for i, month in enumerate(months):

        data.sel(month=month).plot.pcolormesh(ax=axes[i], cmap="coolwarm", 
                                              vmin=vmin, vmax=vmax, extend="both",
                                              add_colorbar=True)

        axes[i].coastlines()
        #axes[i].set_ylabel(month_names[i])

    fig.tight_layout()
    fig.suptitle(title, fontsize=16, y=1.02)

    fig.show()

    # Used suggestions from http://xarray.pydata.org/en/stable/examples/monthly-means.html