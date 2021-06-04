"""
This script pulls out the models, scenarios and runs that contain all the correct
data that we need to run UTCI and saves it in a CSV file.
"""

import numpy as np
import glob
import pandas as pd
from copy import copy
from netCDF4 import Dataset

grid_chunk_size = 40000 * 86 * 365 * 8  # lat-lon x years x days x times per day
# this is the number of data points we'll allow in a model before we break up the run into chunks
# As a guide, SSP scenarios in the BCC model just about complete (51200 lat x lon) whereas CMCC (55296) do not

variables = ['tas', 'huss', 'ps', 'uas', 'vas', 'rsds', 'rsdsdiff', 'rsus', 'rlds', 'rlus']

dirlist = sorted(glob.glob('/badc/cmip6/data/CMIP6/ScenarioMIP/*/*/*/*/3hr/*/g*/latest/'))
insts = [i.split('/')[6] for i in dirlist]
models = [i.split('/')[7] for i in dirlist]
scens = [i.split('/')[8] for i in dirlist]
runs = [i.split('/')[9] for i in dirlist]
vars = [i.split('/')[11] for i in dirlist]

df = pd.DataFrame(list(zip(insts, models, scens, runs, vars)),
               columns = ['institution', 'model', 'scenario', 'run', 'variable'])


dirlist = sorted(glob.glob('/badc/cmip6/data/CMIP6/CMIP/*/*/historical/*/3hr/*/g*/latest/'))
insts = [i.split('/')[6] for i in dirlist]
models = [i.split('/')[7] for i in dirlist]
scens = [i.split('/')[8] for i in dirlist]  # always historical
runs = [i.split('/')[9] for i in dirlist]
vars = [i.split('/')[11] for i in dirlist]

df_hist = pd.DataFrame(list(zip(insts, models, scens, runs, vars)),
               columns = ['institution', 'model', 'scenario', 'run', 'variable'])


models_out = []
scens_out = []
runs_out = []
startyears_out = []
endyears_out = []

for model in df.model.unique():
    for scenario in df[df['model']==model].scenario.unique():
        for run in df[(df['model']==model) & (df['scenario']==scenario)].run.unique():
            vars_exist = df[(df['model']==model) & (df['scenario']==scenario) & (df['run']==run)].variable.tolist()
            count = 0
            for variable in variables:
                if variable in vars_exist and variable in df_hist[(df_hist['model']==model) & (df_hist['run']==run)].variable.tolist():
                    count=count+1
            if count==10:
                datadir = '/badc/cmip6/data/CMIP6/ScenarioMIP/*/%s/%s/%s/3hr/tas/g*/latest/*' % (model, scenario, run)
                filename = glob.glob(datadir)[0]
                nc = Dataset(filename)
                nlat = nc.variables['lat'].shape[0]
                nlon = nc.variables['lon'].shape[0]
                ntime = 86 * 365 * 8
                nchunks = int(np.ceil(ntime * nlat * nlon / grid_chunk_size))
                for ichunk in range(nchunks):
                    startyear = 2015 + ichunk * 86//nchunks
                    endyear = 2015 + (ichunk + 1) * 86//nchunks
                    models_out.append(model)
                    scens_out.append(scenario)
                    runs_out.append(run)
                    startyears_out.append(startyear)
                    endyears_out.append(endyear)
                ntime = 30 * 365 * 8
                nchunks = int(np.ceil(ntime * nlat * nlon / grid_chunk_size))
                for ichunk in range(nchunks):
                    startyear = 1985 + ichunk * 30//nchunks
                    endyear = 1985 + (ichunk + 1) * 30//nchunks
                    models_out.append(model)
                    scens_out.append('historical')
                    runs_out.append(run)
                    startyears_out.append(startyear)
                    endyears_out.append(endyear)


df_out = pd.DataFrame(list(zip(models_out, scens_out, runs_out, startyears_out, endyears_out)), columns=['model', 'scenario', 'run', 'startyear', 'endyear'])
df_out.drop_duplicates(inplace=True, ignore_index=True)

df_out.to_csv('3hr_models.csv')
