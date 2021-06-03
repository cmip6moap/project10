import os
import xarray as xr
from pathlib import Path

data_path = Path("/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg")

def define_path(model, scenario, base_path=data_path, run="r1i1p1f1"):
    ''' Define path to each set of files '''
    path = Path(os.path.join(base_path, model, scenario, run))
    return path

def extract_data(model, scenario, run="r1i1p1f1", base_path=data_path, date_range=None, parameter="utci", resample="D", set_chunks=True):
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

    ## TODO: remove this temporary line
    ##filenames = list(filenames)[0:2]

    if set_chunks:
        data_combined = xr.open_mfdataset(filenames, concat_dim="time", chunks={'time': -1, 'lat': 10, 'lon':10})
        if parameter:
            data_combined = data_combined[parameter]
        if resample:
            data_combined = data_combined.resample(indexer={time_col:resample}).mean()

        print(data_combined)
    else:
        for filename in filenames:

            #filename = next(filenames)

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

    # TODO: At the moment we are reading all the data within a folder and then filtering this afterwards. 
    # Could filter the data files to be read first as this would be more efficient.
    if date_range is not None:
        if len(date_range) == 2:
            data_combined = data_combined.sel(**{time_col:slice(*date_range)})
        else:
            raise ValueError(f"Did not understand input for date_range: {date_range}. Should contain a start and end date only")

    return data_combined

def calc_quantile(ds, quantile=0.95):
    ''' Used to calculate quantile value on a dataset when grouping data '''
    return ds.quantile(quantile, dim="month")

def calc_climatology(data, time_col="time", percentile=None, quantile=0.95, set_chunks=True):
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
    
    Output:
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
        data = data.chunk(chunks={"month":-1})

    if percentile is not None:
        quantile = percentile/100.

    data_quantile = data.groupby("month").map(calc_quantile, quantile = quantile)
    
    time_as_strings = data[time_col].astype('M8[us]').dt.strftime("%Y-%m-%dT%H:%M:%S").sortby(time_col)
    data_quantile.attrs["timeframe"] = f"Climatology over time range: {time_as_strings.values[0]} - {time_as_strings.values[-1]}"
    data_quantile.attrs["quantile"] = f"Quantile calculated: {quantile} (i.e. {quantile*100.:.0f}% percentile)"
    
    return data_quantile

def define_output_name(path, model, scenario, date_range, percentile=95):
    ''' Define output name for netcdf files containing climatology calculated at a specified percentile '''
    filename = os.path.join(path, f"UTCI_climatology_p{percentile:.0f}_{model}_{scenario}_{date_range[0]}-{date_range[1]}.nc")
    return filename

def define_output_name_diff(path, percentile=95):
    ''' Define output name for overall difference utci parameter '''
    filename = os.path.join(path, f"UTCI_climatology_p{percentile:.0f}_difference.nc")
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

def calc_percentile_difference(output_path, write_steps=True, set_chunks=True):
    '''
    Calculate the difference in the 95th percentile across the scenarios for each model.

    This calculate the monthly climatology across two 30 year periods, comparing each ssp 
    scenario to the historical data set.

    historical dataset - from 1985-2014 (based on current data availability)
    ssp scenario - from 2029-2058 (based on current data availability)
    
    models - "BCC-CSM2-MR"
    scenarios - "ssp126", "ssp245", "ssp585"

    Args:
        output_path (str/pathlib.Path):
            Where to write output data from this calculation.
        write_steps (bool) :
            Whether to write out data after the percentile calculation for each model and scenario as well.
            Default = True
    
    Returns:
        xarray.DataArray :
            Difference data for the scenarios wihin each model as one combined DataArray object.
            coords: model, scenario, month, lat, lon

    TODO: Hard-coded at the moment but can add more flexibility in future.
    '''
    # Gathering together details for different models and scenarios we will be comparing
    models = ["BCC-CSM2-MR"]#, "HadGEM3-GC31-LL"]
    ssp_scenarios = ["ssp126", "ssp245", "ssp585"]
    reference_scenario = "historical" # Define folder name for data to use as the reference timeseries

    ref_date_range = ["1985", "2014"]
    scenario_date_range = ["2029", "2058"]

    run = "r1i1p1f1" # Name of lowest level folder (run?)
    percentile=95

    # For each model process the historical data and data for each scenario (assumes same scenarios per model)
    model_diff = []
    for model in models:
        # Calculate 95th percentile for reference (historical data)
        # Check if this already exists first
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
    filename = define_output_name_diff(output_path, percentile=95)
    data_diff.to_dataset().to_netcdf(filename)

    return data_diff

if __name__=="__main__":

    output_path = Path("/home/users/rt17603/acrg_gws/rt17603/CMIP6_output")
    write_steps = True

    data = calc_percentile_difference(output_path, write_steps=write_steps)

    #print(data)

    #average_per_model = data.mean(dim="scenario")
    #average_overall = data.mean(dim=["model","scenario"])
    #print(average_per_model)
    #print(average_overall)
