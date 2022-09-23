## Preparation
```
# Install "snakemake"
mamba create -n snakemake python snakemake
mamba activate snakemake

# Clone this repo
git clone https://github.com/TreesLab/circRNA_features.git
```

### Prepare the reference files
```
cd circRNA_features/resources/
tar zxvf references.tar.gz
```


## Run pipeline
```
snakemake --use-conda --configfile config/human.yaml --resources mem_mb=450000 IO_limit=4 -c 30
```
