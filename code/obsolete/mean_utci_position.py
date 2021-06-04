import argparse
import pandas as pd
from dask.diagnostics import ProgressBar
import xarray as xr

parser = argparse.ArgumentParser(description = "Extract monthly UTCI means for a location")
parser.add_argument('--populated')
parser.add_argument('--row', type=int)
parser.add_argument('--netcdf', nargs="+")
parser.add_argument('--outfile')


args = parser.parse_args()

# For the xarray object
# - select the coordinates for a given city
# - calculate the monthly mean UTCI at those coordinates
# - return a dataframe with utci over time, as well as the city info
def mean_utci_position(utci, pop, row):
    with ProgressBar():
        out = (utci
               .sel(lon=pop['lon'][row], lat=pop['lat'][row], method='nearest')['utci']
               .resample(time="1M")
               .mean()
               .compute()
               .to_dataframe()
               .assign(
                    name=pop['name'][row],
                    adm0name=pop['adm0name'][row],
                    pop=pop['pop'][row],
                    lon_orig=pop['lon'][row],
                    lat_orig=pop['lat'][row]
                ))
    return(out)

populated = pd.read_parquet(args.populated)

l = []
for netcdf in args.netcdf:
    print(netcdf)
    utci =  xr.open_mfdataset(netcdf, parallel=True)
    out = mean_utci_position(utci, populated, args.row)
    out = out.assign(netcdf = netcdf)
    l.append(out)

out = pd.concat(l).reset_index()
out.to_parquet(args.outfile)


