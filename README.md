# LiverDB Web App
Shiny app for exploring the liver RNA-Seq data analyzed by *Ilieva et al. 2022*.

Visit the live app here: [https://rnamedicine.shinyapps.io/liverdb/](https://rnamedicine.shinyapps.io/liverdb/)

See [`preprocess/`](https://github.com/Bishop-Laboratory/LiverDB/tree/main/preprocess) for details on how the data was generated.

# Launching the App Locally
These are the instructions to launch the web app locally via Docker.

Before you get started, make sure:
- You have `git lfs` installed ([see here](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage)) before you clone the repo
    - `liverdb/app_data.rds` should have a file size ~107MB
- You have `docker` installed ([see here](https://docs.docker.com/engine/install/ubuntu/#install-using-the-convenience-script))
    - The Dockerized version of the app was tested with `Docker version 20.10.17, build 100c701`, but any recent-ish version of Docker should work

If you wish to change the host/port number, you'll need to edit lines 22 & 25 of the `Dockerfile`, and adjust the steps below accordingly.

## Steps
1. Build the image. The uncompressed image is about 1.85GB. Building may take a while on the first run

```shell
docker build -t liverdb .
```

2. Spin up a container

```shell
docker run -d --rm -p 3838:3838 liverdb
```

3. Access the web app by entering `http://localhost:3838` in your web browser

4. When you're done, stop the container. You can find the container ID by via the command `docker ps`

```shell
docker stop <CONTAINER_ID>
```