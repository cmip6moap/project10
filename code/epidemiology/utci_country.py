import xarray as xr
import numpy as np
import pandas as pd
from dask.diagnostics import ProgressBar
import argparse

parser = argparse.ArgumentParser(description = "Extract monthly UTCI values for a country")
parser.add_argument('--popfile')
parser.add_argument('--country', type=int)
parser.add_argument('--combinations')
parser.add_argument('--outfile')
args = parser.parse_args()

# args = parser.parse_args(['--country', '125', '--popfile', '../../data/processed_data/populated.parquet', '--combinations', '../../data/processed_data/combinations.parquet', '--outfile', '../../data/processed_data/utci_country_monthly/utci_monthly_125.parquet'])

def country_utci(utci, popinfo, country, model="", scenario="", quantile=""):
    pop = (
        popinfo
        .query('adm0name == "{}"'.format(country))
        .reset_index()
    )
    print(pop.shape[0])
    l = []
    for row in range(pop.shape[0]):
        print(row)
        out = (
            utci
            .sel(lon=pop['lon'][row], lat=pop['lat'][row], method='nearest')['utci']
            .compute()
            .to_dataframe()
            .assign(
                name=pop['name'][row],
                adm0name=country,
                pop=pop['pop'][row],
                lon_orig=pop['lon'][row],
                lat_orig=pop['lat'][row],
                model=model,
                scenario=scenario,
                quantile=quantile
            )
        )
        if isinstance(out.index[0], np.floating):
            out = (out
                .assign(
                    year = [int(str(x)[0:4]) for x in out.index],
                    month = [int(str(x)[4:6].lstrip('0')) for x in out.index]
                )
                .reset_index(drop=True))
        else:
            out = (out
                .assign(
                    year = [x.year for x in out.index],
                    month = [x.month for x in out.index]
                )
                .reset_index(drop=True))
        l.append(out)
    out = pd.concat(l)
    popsum = pop['pop'].sum()
    o = (out
         .assign(utci=out['pop'] * out['utci']/popsum)
         .groupby(['year', 'month', 'adm0name', 'model', 'scenario', 'quantile'])
         .agg({'utci':['sum']})
         .reset_index()
    )
    o.columns = o.columns.droplevel(1)
    return(o)


combinations = pd.read_parquet(args.combinations)
pop = pd.read_parquet(args.popfile)
countries = list(set(pop['adm0name']))
countries.sort()
country = countries[args.country]
print(country)

l = []
for row in range(combinations.shape[0]):
    print(row)
    utci = xr.open_dataset(combinations['datasets'][row])
    l.append(country_utci(utci, pop, country, combinations['model'][row], combinations['scenario'][row], combinations['quantile'][row]))

out = pd.concat(l)
out.to_parquet(args.outfile, compression="gzip")

