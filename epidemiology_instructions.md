# Snakemake pipeline

Author: Gibran Hemani

To setup datasets run:

```
source /gws/pw/j05/cop26_hackathons/bristol/activate-env
snakemake -j 1 -s Snakefile_setup
```

Then to extract UTCI and model mortality run:

```
snakemake -prk \
-j 10 \
--cluster-config slurm.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output}"
```




