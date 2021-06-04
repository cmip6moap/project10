import os
import pandas as pd

os.makedirs("data/processed_data/utci_country_monthly/countries", exist_ok=True)
os.makedirs("job_reports", exist_ok=True)

combinationsfile = "data/processed_data/combinations.parquet.gz"
popfile = "data/processed_data/populated.parquet.gz"
combinations = pd.read_parquet(combinationsfile)
populated = pd.read_parquet("data/processed_data/populated.parquet.gz")
countries = list(set(populated['adm0name']))
countries.sort()
ci = list(range(len(countries)))

rule all:
  input:
    expand("data/processed_data/countries/utci_monthly_{country}.parquet.gz", country=ci)

rule utci_countries:
  input:
    expand("{combinationsfile}", combinationsfile=combinationsfile),
    expand("{popfile}", popfile=popfile)
  output:
    "data/processed_data/countries/utci_monthly_{country}.parquet.gz"
  shell:
    """
      cd code
      python utci_country.py \
        --country {wildcards.country} \
        --popfile ../{popfile} \
        --combinations ../{combinationsfile} \
        --outfile ../{output}
    """

