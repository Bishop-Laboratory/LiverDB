# Preprocess LiverDB data

These are the instructions used to re-build all the datasets used by FibroDB.

Before starting, make sure:
- Your working directory is `preprocess/`
- You have `conda` installed ([see here](https://docs.conda.io/en/latest/miniconda.html))

To generate `app_data.rds` (which houses all the data for the web app), you have two options:
- [Option 1: Generate Processed Data From Raw Files](#option-1-generate-processed-data-from-raw-files)
- [Option 2: Download Processed Data from S3](#option-2-download-processed-data-from-s3)

## Option 1: Generate Processed Data From Raw Files

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
A successful run will generate the directory `raw_counts/` with `.out.tab` files (soft links). For data validation, you can check your generated `.out.tab` files against MD5 hashes in [`metadata/md5sums.txt`](https://github.com/Bishop-Laboratory/LiverDB/blob/main/preprocess/metadata/md5sums.txt)

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
Rscript prepare_app_data.R
```

A successful run will generate `app_data.rds`.

6. Move `app_data.rds` to the app directory

```shell
mv app_data.rds ../liverdb/
```


## Option 2: Download Processed Data from S3

1. Install mamba, create the environment and activate it

```shell
conda install -c conda-forge mamba -y
mamba env create -f environment.yml
conda activate prep
```

2. Download `.csv.gz` files from S3 into the working directory

```shell
wget https://liverdb-data.s3.amazonaws.com/GSE126848_degs.csv.gz
wget https://liverdb-data.s3.amazonaws.com/GSE135251_degs.csv.gz
wget https://liverdb-data.s3.amazonaws.com/GSE126848_gene_exp.csv.gz
wget https://liverdb-data.s3.amazonaws.com/GSE135251_gene_exp.csv.gz
wget https://liverdb-data.s3.amazonaws.com/enrichr_res.csv.gz
```

3. Run `prepare_app_data.R`

```shell
Rscript prepare_app_data.R
```

A successful run will generate `app_data.rds`.

4. Move `app_data.rds` to the app directory

```shell
mv app_data.rds ../liverdb/
```
