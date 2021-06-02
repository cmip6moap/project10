import os
import glob
import xarray as xr
from pathlib import Path
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from .climatology import monthly_cycle

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
scenario = scenarios[0]

# Find all files within relevant folder
full_path = define_path(model, scenario, run=run, 
                        base_path=path)
filenames = full_path.glob("*.nc")

# Define empty variable to use for creating a Dataset
# TODO: See if there is a better way to combine the data to look at a climatology
utci = None

# TODO: Replace this to look at all filenames. Just looking at the first 2 files for now to check this works
# **Temporary loop for now**
for i in range(2):

    filename = next(filenames)

    with xr.open_dataset(filename) as ds:
        utci_one_year = ds["utci"]

    if utci is None:
        utci = utci_one_year
    else:
        utci = xr.concat([utci, utci_one_year], dim="time")

    print(utci)

# May want to make sure the data is explicitly sorted? Not sure which order the generator returns the files.
# May want to filter this to explicitly cover a particular time period (could do this at the filenames stage to avoid reading in unecessary data)

# Use pre-built (old) climatology code to calculate the monthly averages across the whole time period
utci_monthly = monthly_cycle(utci)
utci_monthly.quantile(0.95, dim = ["month", "lat", "lon"])

print(utci_monthly)
