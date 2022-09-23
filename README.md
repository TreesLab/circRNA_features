## Preparation
```
# Install "snakemake"
mamba create -n snakemake python snakemake
mamba activate snakemake

# Clone this repo
git clone https://github.com/TreesLab/circRNA_features.git
cd circRNA_features
```

### Prepare the reference files
```
wget treeslab1.genomics.sinica.edu.tw/circRNA_features/resources.tar.gz
tar zxvf resources.tar.gz
```
Please prepare the reference files described in "resources/README.md".


## Run pipeline
```
snakemake --use-conda --configfile config/human.yaml --resources mem_mb=450000 IO_limit=4 -c 30
```
