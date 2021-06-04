# Project 10: Human heat stress in a warming world

This project looks at how the heat stress metrics are projected to change under different future scenarios. We focus particularly on UTCI.

## Contributors

* Lead: Chris Smith (@chrisroadmap) https://environment.leeds.ac.uk/see/staff/1542/dr-chris-smith
* Charles H. Simpson [@C-H-Simpson](https://github.com/C-H-Simpson)
* Chloe Brimicombe
* Claudia Di Napoli [@cladinapoli](https://github.com/cladinapoli) https://www.reading.ac.uk/geographyandenvironmentalscience/About/Staff/c-dinapoli.aspx
* Gibran Hemani [@explodecomputer](https://github.com/explodecomputer)
* Laila Gohar [@lkgohar](https://github.com/lkgohar)
* Lauren Burton [@eeleb](https://github.com/eeleb)
* Michael Taylor [@patternizer](https://github.com/patternizer)
* Rachel Tunnicliffe [@rt17603](https://github.com/rt17603)
* Robin Lamboll [@rlamboll](https://github.com/rlamboll)
* Seb Cole

## What was done
UTCI in a few scenarios were calculated. Bias correction in this was begun. Various plots of the time evolution of this over the globe were made

### How we approached the problem and why

Traditional extremes of heat and other climate indices defined by the ETCCDI such as monthly maximum temperature or heatwave duration (Zhang et al., 2011) do not directly estimate the impacts on human welfare. Here we compute a metric of human heat stress, the Universal Thermal Climate Index (UTCI; Błazejczyk et al., 2013), using CMIP6 projections. We compare these to wet bulb temperatures across several models. 

### Data we used and how to obtain this

We require 3hrly data including insolation from historical data from 1985, and projections from SSP5-8.5, SSP2-4.5 and SSP1-2.6 to 2100. Most calculations use processed monthly versions of this, focusing on either 95% high values or mean values. Models must include diffuse irradiance to calculate this. Chris collected a few such models and precalculated the UTCI data. Charles Simpson then calculated monthly 95th quantiles and mean data which was also put onto Jasmin. 

### What we did during the hackathon

* Calculate UTCI monthly mean and 95%
* Calculate wet bulb temperatures
* Calculate biases in reproducing historical records
* Investigate the physiological response of bodies to these metrics
* Plot the projections in many and various ways
* Calculate correlations between UTCI and GSAT 

### Outcomes

* Many exciting graphs
* Newfound appreciation for the difficulty of working on Jasmin

## About this repo

There are further `README` files in key directories.

### Key files

* [...]
* [...]
* [...]

### How to reproduce our outputs

1. Find all the data we used on Jasmin (it will probably have moved)
2. Run all of the notebooks

### Repo structure

    .
    ├── notebooks
    │           The Jupyter Notebooks that we created
    │
    ├── code
    │           Any code (Python or otherwise) that we created that doesn't
    │           sit within a Notebook
    │
    ├── results
    │           The key figures that we produced
    │
    ├── data
    │   ├── raw_data
    │   │       Any data we used that didn't come from JASMIN
    │   │
    │   └── processed_data
    │           Any output data that we produced
    │
    ├── environment.yml
    └── environment_frozen.yml
            The libraries and versions that we used

## Next steps for our project

* [...]
* [...]
* [...]
