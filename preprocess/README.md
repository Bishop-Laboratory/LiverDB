# Preprocess LiverDB data

This is the protocol used to re-build all the datasets used by FibroDB.

This protocol assumes:
- Your working directory is `preprocess/`
- You have `conda` installed ([see here](https://docs.conda.io/en/latest/miniconda.html))

## Steps

1. Install mamba, create the environment and activate it

```shell
conda install -c conda-forge mamba -y
mamba env create -f environment.yml
conda activate prep
```

2. Run the Nextflow pipeline

```shell
nextflow run pipeline.nf
```

You might need to edit the constants in `pipeline.nf` depending on how many cores you have available.
A successful run will generate the directory `raw_counts/` with `.out.tab` files (soft links).

3. Run `process_raw_counts.R`

```shell
Rscript process_raw_counts.R
```

A successful run will generate the following files, which are required for the next step:
- `GSE126848_degs.csv.gz`
- `GSE135251_degs.csv.gz`
- `GSE126848_gene_exp.csv.gz`
- `GSE135251_gene_exp.csv.gz`

4. Run `enrichr_res.R`

```shell
Rscript enrichr_res.R
```

A successful run will generate `enrichr_res.csv.gz`.

5. Run `prepare_app_data.R`

```shell
Rscript enrichr_res.R
```

A successful run will generate `app_data.rds`.

6. Move `app_data.rds` to the app directory

```shell
mv app_data.rds ../liverdb/
```

## Raw Counts Validation
You can check the your generated `.out.tab` files against MD5 hashes in [`metadata/md5sums.txt`](https://github.com/Bishop-Laboratory/LiverDB/blob/main/preprocess/metadata/md5sums.txt)
