# Base image https://hub.docker.com/u/rocker/
FROM rocker/shiny:4.2.0

# system libraries of general use
## install debian packages
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev

# copy the renv.lock file for package installation
COPY /liverdb/renv.lock ./renv.lock

# install renv & restore packages
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'renv::consent(provided = TRUE)'
RUN Rscript -e 'renv::restore()'

# copy app files
COPY /liverdb ./app

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]