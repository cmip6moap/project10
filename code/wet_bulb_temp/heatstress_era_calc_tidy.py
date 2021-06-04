"""
Module to calculate Wet Bulb Globe Temperature estimate.
Python code using the Australian Bureau of Meteorolgy's approximation
assuming moderately high radiation level in light wind conditions.
The original formula of Wet Bulb Globe Temperature is:
  WBGT = 0.7 * T_exposed_wet-bulb + 0.2 * T_black_globe + 0.1 * T_screened_dry_bulb
  
The BoM's approximation does not take in to account of solar radiation or wind speed
which would affect the exposed wet bulb temperature and the black globe temperature
in reality and hence the value of WBGT. 
  WBGT_approx = 0.567 * T_dry_bulb + 0.393 * vapour_pressure + 3.94

Taken from   http://www.bom.gov.au/info/thermal_stress/

These values have been corrected for bias found by Damian Wilson and colleagues.
However the above estimate is for daily maximum temperatures - so should take the maximum for each day.

This will calculate WBGT using ERA-5 data and requires, temperature ('2t'), dewpoint ('2d') 

Adapted from Chris Smith's excellent scripts for UTCI

Written by Laila Gohar 04062021 laila.gohar@metoffice.gov.uk

-----------------------------------------------
"""

import numpy as np
from satvap_calc_array_era import GetSatVapPres
from satvap_calc_array_era import GetTKelvinFromTCelsius
#from humidity import calc_saturation_vapour_pressure, specific_to_relative
EPSILON = 0.62198  # ratio of water vapour to dry air molecular weights


def calc_vapour_pressure(tas,rel_hum):
    """Calculate vapour pressure 
       vapour pressure - e
       relative humidity - rh
       dry bulb temperature - tas
       e = rh / 100 × 6.105 × exp ( 17.27 × Ta / ( 237.7 + Ta ) )

    """
    rh = rel_hum
    ta = tas - 273.15
    partexp = 17.27*ta/(237.7 + ta)

    result = (rh/100)*6.105*np.exp(partexp)

    return result

# TODO:
# - throw warning if any of the input parameters are out of range the relationships
#   were designed for


def wet_bulb_globe_temperature_era_BoM(base):

    """Calculate Wet Bulb Globe Temperature approximations.
    Parameters
    ----------
        base : dict of array_like
            dict containing CMIP-style variables, which should contain the following
            keys:
            2t     : near-surface air temperature, K
            2d    : dew point temperature, K
   Returns
    -------
        wbgt : array_like
            Wet Bulb Globe Temperature value, K
    Raises
    ------
        ValueError:
            if input variables do not match what is required by the model.
    """
    # this appears in APRP code twice: refactor target
    #check_vars = ["tas"]
    check_vars = ["2t", "2d"]
    for check_var in check_vars:
        if check_var not in base.keys():
            raise ValueError("%s not present in %s" % (check_var, "base"))



    # set default pressure if not in input
    if "ps" not in base.keys():
        base["msl"] = 101325
    # allow list input: convert to array
    base["2t"] = np.asarray(base["2t"])
    
    # turn off pylint warnings on short names, because the calculation gets stupid
    ta = base["2t"] - 273.15  # pylint: disable=invalid-name

    base["2d"] = np.asarray(base["2d"])
    vapour_pressure = GetSatVapPres(base["2d"])
    # Vapour pressure from this formula is in pascals. need hectorpascals for wBGT BoM formula
    #saturation_vapour_pressure = calc_saturation_vapour_pressure(base["tas"])


    # WBGT approximation:
    wbgt_bom = 0.567 * ta + 0.393 * vapour_pressure*1e-2 + 3.94

    #**** Correct for biases in the approximation as derived by Damian Wilson and colleagues
    #***************************************************************************************
    #Bias for condition 1: WBGT < 24.82degC 
    wbgt_dum = wbgt_bom.copy()
    wbgt_bom1 = np.where(wbgt_dum < 24.82, wbgt_dum-1.32,wbgt_dum) 
    
    #Bias for condition 2 : 24.82degC < WBGT < 40.0degC
    bias = 0.013497 * wbgt_bom1**2 - 0.6700 * wbgt_bom1 + 9.64
    wbgt_bom2 =  np.where((wbgt_bom1 > 24.82) & (wbgt_bom1 < 40.), wbgt_bom1-bias,wbgt_bom1) 

    #Bias for condition 3: 40.0degC < WBGT
    bias1 = 0.4098 * wbgt_bom2 - 11.96
    result = np.where(wbgt_bom2 >= 40.0,wbgt_bom2-bias1, wbgt_bom2) 

    return result
