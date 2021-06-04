import os
import pandas as pd

os.makedirs("data/processed_data/utci_country_monthly/france", exist_ok=True)
os.makedirs("data/processed_data/utci_country_monthly/all", exist_ok=True)
os.makedirs("job_reports", exist_ok=True)


models = ['HadGEM3-GC31-LL', 'BCC-CSM2-MR']
runs = ['historical', 'ssp126', 'ssp245', 'ssp585']
projectdir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg'

utci_dirs = [
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/historical/r1i1p1f3/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/ssp126/r1i1p1f3/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/ssp245/r1i1p1f3/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/HadGEM3-GC31-LL/ssp585/r1i1p1f3/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/BCC-CSM2-MR/historical/r1i1p1f1/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/BCC-CSM2-MR/ssp126/r1i1p1f1/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/BCC-CSM2-MR/ssp245/r1i1p1f1/*.nc',
	'/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/BCC-CSM2-MR/ssp585/r1i1p1f1/*.nc'
]

models = '"' + '" "'.join(utci_dirs) + '"'

populated = pd.read_parquet("data/processed_data/populated.parquet.gz")
populated

france = pd.read_parquet("data/processed_data/france.parquet.gz")
france

franceindex = list(range(france.shape[0]))
allindex = list(range(populated.shape[0]))

rule all:
	input:
		expand("data/processed_data/utci_country_monthly/france/{franceindex}.parquet.gz", franceindex=franceindex),
		expand("data/processed_data/utci_country_monthly/all/{allindex}.parquet.gz", allindex=allindex)


rule utci_france:
	input:
		"data/processed_data/france.parquet.gz"
	output:
		"data/processed_data/utci_country_monthly/france/{franceindex}.parquet.gz"
	conda:
		"jaspy3.7-m3-4.9.2-r20210320"
	shell:
		"""
		cd code
		python mean_utci_position.py \
		--netcdf {models} \
		--row {wildcards.franceindex} \
		--populated ../data/processed_data/france.parquet.gz \
		--outfile ../{output}
		"""

rule aggregate_france:
	input:
		expand("data/processed_data/utci_country_monthly/france/{franceindex}.parquet.gz", franceindex=franceindex)
	output:
		"france_done"
	shell:
		"touch {output}"

rule utci_all:
	input:
		"data/processed_data/populated.parquet.gz",
		rules.aggregate_france.output
	output:
		"data/processed_data/utci_country_monthly/all/{allindex}.parquet.gz"
	conda:
		"jaspy3.7-m3-4.9.2-r20210320"
	shell:
		"""
		cd code
		python mean_utci_position.py \
		--netcdf {models} \
		--row {wildcards.allindex} \
		--populated ../data/processed_data/populated.parquet.gz \
		--outfile ../{output}
		"""

	


