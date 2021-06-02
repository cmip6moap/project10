import os
import glob
import xarray as xr
from pathlib import Path
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from climatology import monthly_cycle

def define_path(model, scenario, base_path=Path("/home/users/rt17603/cop26_hackathons/project10/utci_projections_1deg"),
                run="r1i1p1f1"):
    ''' Define path to each set of files '''
    path = Path(os.path.join(base_path, model, scenario, run))
    return path

path = Path("/home/users/rt17603/cop26_hackathons/project10/utci_projections_1deg")

# Gathering together details for different models and scenarios we will be comparing
models = ["BCC-CSM2-MR", "HadGEM3-GC31-LL"]
scenarios = ["ssp126", "ssp245", "ssp585"]
reference_scenario = "historical" # Define folder name for current data to use as the reference timeseries

run = "r1i1p1f1" # Name of lowest level folder

# TODO: Wrap this into loop etc. to run over different models and scenarios
# **Temporary lines for now**
model = models[0]
scenario = reference_scenario

# Find all files within relevant folder
full_path = define_path(model, scenario, run=run, 
                        base_path=path)
filenames = full_path.glob("*.nc")

# Define empty variable to use for creating a Dataset
# TODO: See if there is a better way to combine the data to look at a climatology
utci = None

# TODO: Replace this to look at all filenames. Just looking at the first 2 files for now to check this works
# **Temporary loop for now**
# Could replace with open_mfdatasets(filenames, concat_dim="time") but then would be slower if we'll continue to apply the resample
# - could see if preprocess input would allow you to apply e.g. a resample operation prior to concatenation
for i in range(2):

    filename = next(filenames)

    with xr.open_dataset(filename) as ds:
        utci_daily = ds["utci"].resample(indexer={"time":"D"}).mean() # May want to remove this but is reducing data quantity for now

    if utci is None:
        utci = utci_daily
    else:
        utci = xr.concat([utci, utci_daily], dim="time")

    print(utci)

date_range = ["1985", "2014"]

utci = utci.sortby("time")
utci = utci.sel(time=slice(*date_range))

# May want to make sure the data is explicitly sorted? Not sure which order the generator returns the files.
# May want to filter this to explicitly cover a particular time period (could do this at the filenames stage to avoid reading in unecessary data)

# Use pre-built (old) climatology code to calculate the monthly averages across the whole time period
utci_ref_95 = monthly_cycle(utci, quantile=0.95)

print(utci_ref_95)
