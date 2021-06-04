"""
Developed by: Rachel Tunnicliffe

This module is used to read and process 3-hourly UTCI output derived from CMIP6 models.

Function call
=============

The main wrapper function is::

    calc_percentile_difference()

This calculates the percentile values (default 95th) for each scenario within a given 
time period (default 30 years from 2071-2100) and subtracts the same percentile from
historical data (default 30 years from 1985-2014). It does this monthly over the 30 year
time period and produces a DataArray (output netcdf file) containing 12 monthly values for
each lat, lon, model and scenario.

This can be called as::

    calc_percentile_difference(output_path)

where `output_path` specifies where to write the output data file.

This file can be run as::

    $ python utci_diff.py

This will call the __main__ block below and output to the output_path specified there
(to our project10 gws for the hackathon).

Note: ** At the moment this is hardwired to expect certain models and scenarios but this
can be updated either within the function or additional inputs added to allow this to
be updated **

Expected data structure and files
=================================

Expected directory structure for these files::

    model/scenario/run_name/

e.g. BCC-CSM2-MR/historical/r1i1p1f1/

Example filename::

    utci_3hr_BCC-CSM2-MR_historical_r1i1p1f1_gn_198501010300-198601010000.nc

The date string of the form YYYYMMDDhhmmss-YYYYMMDDhhmmss (e.g. 198501010300-198601010000)
is expected and used when extracting the files with an input date range. This is also
expected as a netcdf (.nc) file.

File should contain (at least) the "utci" variable with (time, lat, lon) dimensions.
"""


import os
import re
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

data_path = Path("/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg")

def define_path(model, scenario, base_path=data_path, run="r1i1p1f1"):
    ''' Define path to each set of 3hr UTCI input files '''
    path = Path(os.path.join(base_path, model, scenario, run))
    return path

def extract_files(filenames, start, end):
    '''
    Extract filenames within a given date range.

    Expect filename strings of the format:
        *YYYYMMDDhhmmss-YYYYMMDDhhmmss*.nc
        e.g. ./utci_3hr_BCC-CSM2-MR_historical_r1i1p1f1_gn_198501010300-198601010000.nc

    Input:
        filenames (list) :
            List of filenames extracted from a folder
        start, end (str) :
            Start and end date to use for filtering the file names.
            Will definitely work when start and end is included as "YYYY"
            Should also work for other recognised pandas formats e.g.
            "YYYY-MM-DD".
    
    Returns:
        list :
            List of filenames within the date range specified. Filtered
            from input filenames.

        ValueError:
            No files are found in that date range

    TODO: Only start date extracted from file name is used at present to
    filter the filenamea. Could also use end date.
    '''

    # dateformat example "198501010300-198601010000"

    # Expect input start and end date as string e.g. "2013" for now
    start = pd.to_datetime(start)#, format="%Y")
    end = pd.to_datetime(end)#, format="%Y")

    filenames_match = []
    for filename in filenames:
        filename = str(filename)
        try:
            re_str = "\d{12}[-]\d{12}"
            d = re.search(re_str, filename)
            d = d.group() # Extract value from regular expression compiler
        except AttributeError:
                pass
        else:
            s = d.split('-')[0]
            s = pd.to_datetime(s, format="%Y%m%d%H%M%S")
            # TODO: Could also incorporate a check for end as well as start but should be ok with this data
            # e = d.split('-')[1]
            # e = pd.to_datetime(e, format="%y%m%d%H%M%S")
            if (s >= start and s <= end):
                filenames_match.append(filename)

    if not filenames:
        raise ValueError(f"No filenames found for date range: {start} - {end}")

    return filenames_match

def extract_data(model, scenario, run, base_path=data_path, date_range=None, parameter="utci", resample="D", set_chunks=False):
    '''
    Extract data files for the model and scenario (based on expected directory structure):

        model/scenario/run/*.nc

    Args:
        model (str) :
            Name of model (folder name) e.g. BCC-CSM2-MR
        scenario (str) :
            Name of the scenario (folder name) e.g. ssp126
        run (str, optional) :
            Run name(?) (folder name) e.g. r1i1p1f1
            Default = "r1i1p1f1"
        base_path (str, optional) :
            Base level path for directory structure containing the input files
        date_range (list/None, optional) :
            Start and end date to use to filter the data. At the moment this 
            expects a two-item list containing strings or datetime objects to
            use to slice the time dimension on the input data.
            Default = None
        parameter (str/None, optional) :
            Extract one parameter contained within the Dataset. Set to None to extract
            all parameters.
            Default = "utci"
        resample (str/None, optional) :
            Resample the data along the time axis to reduce the data frequency (mean).
            This needs to match pandas/xarray frequency aliases.
            Set to None to not resample the data.
            Default = "D" (daily)

    Returns:
        xarray.DataArray / xarray.DataSet :
            Extracted and concatanated data for the model / scenario inputted.
            This will be an xarray.DataArray if parameter has been input; xarray.Dataset otherwise
    '''
    
    full_path = define_path(model, scenario, run=run, 
                            base_path=base_path)
    filenames = full_path.glob("*.nc")

    time_col = "time"

    if date_range is not None:
        filenames = extract_files(filenames, date_range[0], date_range[1])

    if set_chunks:
        data_combined = xr.open_mfdataset(filenames, concat_dim="time", chunks={'time': -1, 'lat': 10, 'lon':10})
        if parameter:
            data_combined = data_combined[parameter]
        if resample:
            data_combined = data_combined.resample(indexer={time_col:resample}).mean()
        
        print("Reading data complete.")
    else:
        data_combined = None
        for filename in filenames:

            with xr.open_dataset(filename) as ds:
                if parameter and resample:
                    data_input = ds[parameter].resample(indexer={time_col:resample}).mean() # May want to remove this but is reducing data quantity for now
                elif parameter:
                    data_input = ds[parameter]
                else:
                    data_input = ds

            if data_combined is None:
                data_combined = data_input
            else:
                data_combined = xr.concat([data_combined, data_input], dim=time_col)

    data_combined = data_combined.sortby(time_col)

    if date_range is not None:
        if len(date_range) == 2:
            # Even though input files have been filtered based on the date range specified, 
            # make sure we explicitly filter the data for the dates we want
            data_combined = data_combined.sel(**{time_col:slice(*date_range)})
        else:
            raise ValueError(f"Did not understand input for date_range: {date_range}. Should contain a start and end date only")

    return data_combined

def calc_quantile(ds, quantile=0.95):
    ''' Used to calculate quantile value on a dataset when grouping data '''
    return ds.quantile(quantile, dim="month")

def calc_climatology(data, time_col="time", percentile=None, quantile=0.95, set_chunks=False):
    '''
    Calculate the monthly cycle across time period of input object. Find the quantile for each month.
    
    Args:
        data (xarray.Dataset / xarray.DataArray) :
            Data object to look at seasonal cycle. Must contain a time coordinate.
        time_col (str) :
            Name of time co-ordinate in dataset.
            Default = "time".
        percentile (int/float) :
            Percentile to apply for each month (0-100). If this is specified this will superced the quartile value supplied.
            Default = None (use default quantile value instead)
        quantile (float) :
            Quantile value to apply for each month (0.0-1.0). This is superceded by percentile if specified.
            Default = 0.95
        set_chunks (bool, optional):
            Set to True if chunks are present in the input Dataset already. This will then be applied
            to the new dimension created along the time axis.
            Not aure what this will do otherwise!
            Default = False

    Returns:
        xarray.Dataset / xarray.DataArray
            Data with the quantile calculated for each repeated month.
            "time" dimension is replaced by "month".
            "time" coordinate still contained within Dataset / DataArray but no longer attached to any data variables.
    '''
    months = data[time_col].dt.month

    ## TODO: OLD METHOD BASED ON PREVIOUS CODE - CAN PROBABLY BE IMPROVED
    data = data.assign_coords(**{"month":(time_col, months)})
    data = data.swap_dims({time_col:"month"})
    
    if set_chunks:
        # Set chunks along the month axis expliticly
        data = data.chunk(chunks={"month":-1})

    if percentile is not None:
        quantile = percentile/100.

    data_quantile = data.groupby("month").map(calc_quantile, quantile = quantile)
    
    # Extracting the date as a string was being awkward around different datetime formats for the different models
    # Including this try statement allowed date labels to be extracted from both BCC-CSM2 and HadGEM3.
    # Note: Can remove these lines if they're not working as they're only extracting the dates
    # to save to the attributes.
    try:
        time_as_strings = data[time_col].dt.strftime("%Y-%m-%dT%H:%M:%S").sortby(time_col)
    except ValueError:
        time_as_strings = data[time_col].astype('M8[us]').dt.strftime("%Y-%m-%dT%H:%M:%S").sortby(time_col)
    data_quantile.attrs["timeframe"] = f"Climatology over time range: {time_as_strings.values[0]} - {time_as_strings.values[-1]}"
    data_quantile.attrs["quantile"] = f"Quantile calculated: {quantile} (i.e. {quantile*100.:.0f}% percentile)"
    
    return data_quantile

def define_output_name(path, model, scenario, date_range, percentile=95):
    ''' Define output name for netcdf files containing climatology calculated at a specified percentile '''
    filename = os.path.join(path, f"UTCI_climatology_p{percentile:.0f}_{model}_{scenario}_{date_range[0]}-{date_range[1]}.nc")
    return filename

def define_output_name_diff(path, date_range1, date_range2, percentile=95):
    ''' Define output name for overall difference in utci parameter '''
    filename = os.path.join(path, f"UTCI_climatology_p{percentile:.0f}_difference_{date_range1[0]}-{date_range1[1]}_{date_range2[0]}-{date_range2[1]}.nc")
    return filename

def read_previous_calculation(path, model, scenario, date_range, parameter = "utci", percentile = 95):
    '''
    Check if a file already exists for the given set up and read this data if present.

    Args:
        path (str) :
            Output path where calculation files have been / will be written
        model (str) :
            Name of model e.g. BCC-CSM2-MR
        scenario (str) :
            Name of the scenario e.g. ssp126
        date_range (list):
            Start and end date to use to filter the data. At the moment this 
            expects a two-item list containing strings or datetime objects to
            use to slice the time dimension on the input data.
        parameter (str, optional) :
            Parameter name contained within the Dataset to be extracted.
            Default = "utci"
        percentile (int/float, optional) :
            Percentile value which was used in previous calculation
            Default = 95
    
    Returns:
        xarray.DataArray / None :
            If data is present, Dataset will be read and parameter value extracted. Returned as a DataArray
            Otherwise this will return None
    '''
    out_filename = define_output_name(path, model, scenario, date_range, percentile=percentile)
    if os.path.exists(out_filename):
        with xr.open_dataset(out_filename) as ds:
            if parameter:
                data = ds[parameter]
            else:
                data = ds
        return data
    else:
        return None

def calc_percentile_difference(output_path, ref_date_range = ["1985", "2015"], scenario_date_range=["2071", "2101"], write_steps=True, set_chunks=False):
    '''
    Calculate the difference in the 95th percentile across the scenarios for each model.

    This calculate the monthly climatology across two periods (set by ref_date_range and scenario_date_range),
    comparing each ssp scenario to the historical data set.

    models - "BCC-CSM2-MR", "HadGEM3-GC31-LL"
    scenarios - "ssp126", "ssp245", "ssp585"

    Args:
        output_path (str/pathlib.Path):
            Where to write output data from this calculation.
        ref_date_range (list, optional) :
            Start and end date (end non-inclusive) to extract and use for the reference dataset.
            Should be specified in a recognised pandas date format e.g. "YYYY" or "YYYY-MM-DD"
            Default = ["1985", "2015"] (1985-2014 years included)
        scenario_date_range (list, optional) :
            Start and end date (end non-inclusive) to extract and use for the scenario datasets.
            Should be specified in a recognised pandas date format e.g. "YYYY" or "YYYY-MM-DD"
            Default = ["2071", "2100"] (2071-2100 years included)
        write_steps (bool, optional) :
            Whether to write out data after the percentile calculation for each model and scenario as well.
            Default = True
        set_chunks (bool, optional) :
            Whether to explicitly specify chunks for dask to use. If set to True this will set
            chunks = {"time": -1, "lat": 10, "lon": 10} which means this is not chunked along
            the time dimension. This is preferable when performing operations along the time
            dimension.
            However, this did cause some issues when running on JASMIN sci nodes and made this very slow
            / used a lot of memory (dask does sometimes get memory leaks). So see what works.
            Default = False
    
    Returns:
        xarray.DataArray :
            Difference data for the scenarios wihin each model as one combined DataArray object.
            coords: model, scenario, month, lat, lon

    TODO: Some models and scenarios are hardcoded at the moment but can add these as arguments to allow for
    more flexibility if needed.
    '''

    # Gathering together details for different models and scenarios we will be comparing
    models = ["BCC-CSM2-MR", "HadGEM3-GC31-LL"]
    ssp_scenarios = ["ssp126", "ssp245", "ssp585"]
    reference_scenario = "historical" # Define folder name for data to use as the reference timeseries
    runs = {"BCC-CSM2-MR":"r1i1p1f1", "HadGEM3-GC31-LL":"r1i1p1f3"} # Name of lowest level folder (run?)

    #ref_date_range = ["1985", "2014"]
    #scenario_date_range = ["2029", "2058"]
    
    percentile=95

    # For each model process the historical data and data for each scenario (assumes same scenarios per model)
    model_diff = []
    for model in models:
        # Calculate 95th percentile for reference (historical data)
        # Check if this already exists first
        run = runs[model]
        utci_ref_95 = read_previous_calculation(output_path, model, reference_scenario, ref_date_range, parameter="utci")
        
        if utci_ref_95 is None:
            utci_ref = extract_data(model, reference_scenario, run=run, date_range=ref_date_range, parameter="utci", resample="D",
                                    set_chunks=set_chunks)
            utci_ref_95 = calc_climatology(utci_ref, percentile=percentile, set_chunks=set_chunks)

            if write_steps:
                out_filename = define_output_name(output_path, model, reference_scenario, ref_date_range)
                utci_ref_95.to_dataset().to_netcdf(out_filename)

        # Calculate 95th percentile for each ssp scenario
        # Check if this already exists first
        scenario_diff = []
        for scenario in ssp_scenarios:

            utci_scenario_95 = read_previous_calculation(output_path, model, scenario, scenario_date_range, parameter="utci")

            if utci_scenario_95 is None:
                utci_scenario = extract_data(model, scenario, run=run, date_range=scenario_date_range, parameter="utci", resample="D",
                                             set_chunks=set_chunks)
                utci_scenario_95 = calc_climatology(utci_scenario, percentile=percentile, set_chunks=set_chunks)

                if write_steps:
                    out_filename = define_output_name(output_path, model, scenario, scenario_date_range)
                    utci_scenario_95.to_dataset().to_netcdf(out_filename)
            
            # Calculate the difference between each scenario and the reference and collect together
            utci_diff_95 = utci_scenario_95 - utci_ref_95
            scenario_diff.append(utci_diff_95)
        
        # Concatenate these together into one DataArray with scenario as a dimension (and coordinate)
        scenario_data = xr.concat(scenario_diff, dim="scenario")
        scenario_data = scenario_data.assign_coords({"scenario":ssp_scenarios})
        model_diff.append(scenario_data)
    
    # Concatenate these together again into one DataArray with model as a dimension (and coordinate)
    data_diff = xr.concat(model_diff, dim="model")
    data_diff = data_diff.assign_coords({"model":models})

    # TODO: May want to add this line to make this data output more CF-compliant - dimension order (T, Y, X, ...)
    # data_diff = data_diff.transpose(["month", "lat", "lon", "model", "scenario"])
    
    # Write difference DataArray out to file
    filename = define_output_name_diff(output_path, ref_date_range, scenario_date_range, percentile=95)
    data_diff.to_dataset().to_netcdf(filename)

    return data_diff

if __name__=="__main__":

    output_path = Path("/gws/pw/j05/cop26_hackathons/bristol/project10/percentile_diff_outputs")

    data = calc_percentile_difference(output_path)
