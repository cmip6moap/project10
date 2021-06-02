#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:28:40 2019

@author: rt17603
"""

import numpy as np
import pandas as pd
import xarray as xr
from collections import OrderedDict
import pdb


def mean_month(ds,time_col="time",):
    '''
    Calculate mean across month axis only for each data variable. Do not collapse other axes.
    Used on Grouper object created by ds.groupby("month").
    
    Args:
        ds (xarray.Dataset)
            Emissions data 
        time_col (str)
            Name of the time coordinate
            Defaults to 'time'
    
    Returns:
        xarray.Dataset
            Average emissions per month
            "time" dimension is replaced by integer "month" coordinates.
    '''
    
    if isinstance(ds, xr.core.dataset.Dataset):
        
        data_vars = list(ds.data_vars)
        
        for dv in data_vars:
            if 'month' in ds[dv].dims:
                    ds[dv] = ds[dv].mean(dim="month",keep_attrs=True)
        ds[time_col] = ds[time_col].swap_dims({"month":time_col})
            
        ds["month"]  = ds["month"].mean(dim="month", dtype=np.int)

    else:
        if "month" not in ds.dims:
            # If only one value is grouped, no "month" dims will be present but may still exist within the 
            # DataArray.
            # In this case don't need to average but do want to make sure "month" was a part of this 
            # DataArray before grouping. Expect KeyError to be raised otherwise.
            # Seems like an xarray bug as should really be labelled as a dim.
            check = ds["month"]
        else:
            ds = ds.mean(dim="month",keep_attrs=True)

    return ds

def calc_quantile(ds, quantile=0.95):
    '''
    '''
    return ds.quantile(quantile, dim="month")


def monthly_cycle(ds, time_col="time", quantile=0.95):
    '''
    Calculate the monthly cycle across time period of input object.
    
    Args:
        ds (xarray.Dataset)
            Dataset object to look at seasonal cycle. Must contain a time coordinate.
        time_col (str) 
            Name of time co-ordinate in dataset.
            Default = "time".
    
    Output:
        xarray.Dataset
            Dataset with each variable averaged for each repeated month.
            "time" dimension is replaced by "month".
            "time" coordinate still contained within dataset but no longer attached to any data variables.
    '''
    time = ds[time_col].values
    months = ds[time_col].dt.month

    # if isinstance(ds, xr.core.dataset.Dataset):
    #     dims_list = []
    #     for dv in ds.data_vars:
    #         dv_dims = list(ds[dv].dims)
    #         if time_col in dv_dims:
    #             i = dv_dims.index(time_col)
    #             dv_dims[i] = "month"
                
    #         dims_list.append(dv_dims)
    # else:
    #     dims_list = [list(ds.dims)]
    #     if time_col in ds.dims:
    #         i = ds.dims.index(time_col)
    #         dims_list[0][i] = "month"

    ds_new  = ds.copy(deep=True)
    
    ds_new = ds_new.assign_coords(**{"month":(time_col, months)})
    ds_new = ds_new.swap_dims({time_col:"month"})
    
    ds_quantile = ds_new.groupby("month").map(calc_quantile, quantile = quantile)
    
    #ds_mean = group.map(mean_month,**{"time_col":time_col})
    
    # if isinstance(ds,xr.core.dataset.Dataset):
    #     for i,dv in enumerate(ds.data_vars):
    #         if 'month' in ds_mean[dv].dims and 'month' not in dims_list[i]:
    #             dims_list[i].insert(1, 'month')
    #         ds_mean[dv] = ds_mean[dv].transpose(*dims_list[i])
    # else:
    #     ds_mean = ds_mean.transpose(*dims_list[0])
    

    time_as_strings = ds[time_col].dt.strftime("%Y-%m-%dT%H:%M:%S").sortby(time_col)
    ds_quantile.attrs["timeframe"] = f"Climatology over time range: {time_as_strings.values[0]} - {time_as_strings.values[-1]}"
    ds_quantile.attrs["quantile"] = f"Quantile calculated: {quantile} (i.e. {quantile*100.:.0f}% percentile)"
    
    return ds_quantile
