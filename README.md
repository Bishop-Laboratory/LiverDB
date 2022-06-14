# LiverDB Web App
Shiny app for exploring the liver RNA-Seq data analyzed by *Ilieva et al. 2022*.

Visit the live app here: [https://rnamedicine.shinyapps.io/liverdb/](https://rnamedicine.shinyapps.io/liverdb/)

See [`preprocess/`](https://github.com/Bishop-Laboratory/LiverDB/tree/main/preprocess) for details on how the data was generated.

## Launch the App Locally

These are the instructions used to locally launch the LiverDB web app.

Before starting, make sure:
- You have generated `app_data.rds` and moved it to `liverdb/` (see [`preprocess/`](https://github.com/Bishop-Laboratory/LiverDB/tree/main/preprocess))

### Steps

1. Install mamba, create the environment, activate it and 

```shell
conda install -c conda-forge mamba -y
mamba env create -f liverdb.yml
conda activate liverdb
R -e "install.packages('prompter', repos = 'http://cran.us.r-project.org')"
```

2. Launch the app

```shell
Rscript runApp.R 5566  # Port number
```
