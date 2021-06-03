import os
import pandas as pd

os.makedirs("data/processed_data/utci_country_monthly/france", exist_ok=True)
os.makedirs("data/processed_data/utci_country_monthly/all", exist_ok=True)
os.makedirs("job_reports", exist_ok=True)


models = ['HadGEM3-GC31-LL', 'BCC-CSM2-MR']
runs = ['historical', 'ssp126', 'ssp245', 'ssp585']
projectdir = '/gws/pw/j05/cop26_hackathons/bristol/project10/utci_projections_1deg/'

populated = pd.read_parquet("data/processed_data/populated.parquet.gz")
populated

france = pd.read_parquet("data/processed_data/france.parquet.gz")
france

franceindex = list(range(france.shape[0]))
allindex = list(range(populated.shape[0]))

rule all:
	input:
		expand("data/processed_data/utci_country_monthly/france/{franceindex}_{model}_{run}.parquet.gz", franceindex=franceindex, model=models, run=runs),
		expand("data/processed_data/utci_country_monthly/all/{allindex}_{model}_{run}.parquet.gz", allindex=allindex, model=models, run=runs)


rule utci_france:
	input:
		"data/processed_data/france.parquet.gz"
	output:
		"data/processed_data/utci_country_monthly/france/{franceindex}_{model}_{run}.parquet.gz"
	conda:
		"jaspy3.7-m3-4.9.2-r20210320"
	shell:
		"""
		cd code
		python mean_utci_position.py \
		--netcdf "{projectdir}/{wildcards.model}/{wildcards.run}/r1i1p1f*/*nc" \
		--row {wildcards.franceindex} \
		--populated ../data/processed_data/france.parquet.gz \
		--outfile ../{output}
		"""

rule utci_all:
	input:
		"data/processed_data/populated.parquet.gz"
	output:
		"data/processed_data/utci_country_monthly/all/{allindex}_{model}_{run}.parquet.gz"
	conda:
		"jaspy3.7-m3-4.9.2-r20210320"
	shell:
		"""
		cd code
		python mean_utci_position.py \
		--netcdf "{projectdir}/{wildcards.model}/{wildcards.run}/r1i1p1f*/*nc" \
		--row {wildcards.allindex} \
		--populated ../data/processed_data/populated.parquet.gz \
		--outfile ../data/{output}
		"""

	


