import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def plot_diff_maps(data, title="UTCI difference between ..."):

    fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(14,12), subplot_kw={"projection":ccrs.Robinson()})
    axes = axes.flatten()
    print(axes.shape)

    months = data.month

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    for i, month in enumerate(months):
        print(i, month)
        data.sel(month=month).plot.pcolormesh(ax=axes[i], cmap="twilight_shifted")

        axes[i].set_ylabel(month_names[i])

    fig.tight_layout()
    fig.suptitle(title, fontsize=16, y=1.02)

    # Adapted suggestions from http://xarray.pydata.org/en/stable/examples/monthly-means.html