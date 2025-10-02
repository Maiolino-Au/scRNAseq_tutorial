# Single Cell RNA Sequencing (scRNAseq) with Seurat

Maiolino Aurelio - September 2025

This is still a work in progress

# Index

0. [Introduction](#introduction---work-in-progress)
1. [Docker: the Working Environment](#docker-the-working-environment)
    1. [The basis: how to build and use a Docker image](#the-basis-how-to-build-and-use-a-docker-image)
        1. [The Dockerfile](#the-dockerfile)
        2. [Building the Docker image](#building-the-docker-image)
        3. [Running a Docker container](#running-a-docker-container)
        4. [Saving changes to an image](#saving-changes-to-an-image)
        5. [Usefull Docker commands](#usefull-docker-commands)
    2. [Preexisting Docker images](#preexisting-docker-images)
    3. [Credo: a better and faster way to build Docker images](#credo-a-better-and-faster-way-to-build-docker-images---work-in-progress)
2. [The basis](#the-basis)
    1. [The basis](#the-basis)

 
# Introduction - WORK IN PROGRESS




\newpage

# Docker: the Working Environment
To ensure a consistene and reproducible environement for the analysis, it is reccomended to use a Docker container. Docker allows to create a virtual environemnt that contains all the necessary software and dependencies for the analysis. This way, you can avoid issues related to different versions of R, Seurat, and other packages. With proper version control, you can ensure that the analysis will run smoothly and produce consistent results, even years after the initial analysis.

Creating a Docker means buildng what is called a Docker image. This is formed by various layers, each representing a step in the installation process. The instructions for building the image are written in a file called `Dockerfile` (without any file extension). This text file contains a series of commands that specify the base image, install necessary packages, and set up the environment. The Dockerfile is used to create a Docker image, which can then be run as a container. The container is an isolated environment that contains all the necessary software and dependencies for the analysis. Everythign that happens inside the coontainer will remain inside the single instance and will not affect the imgage from which it was created, unlsess you explicitly save - commit - the changes to a new image.

## The basis: how to build and use a Docker image

Before starting to use docker it is necessary to install it on your system. Docker is available for Windows, macOS, and Linux. The installation process varies slightly depending on the operating system:

- On windows, Docker is installed as a desktop application. It allows you to run containers on your local machine. The Docker Desktop application provides a user interface to manage your containers and images, making it easier to work with Docker. Docker Desktop can be downloaded from the [Docker website](https://www.docker.com/products/docker-desktop/).

- On macOS, Docker is also installed as a desktop application. It provides a similar user interface as on Windows, allowing you to manage your containers and images. Docker Desktop for Mac can be downloaded from the [Docker website](https://www.docker.com/products/docker-desktop/).

- On Linux, Docker is installed as a command-line tool. It allows you to run containers and manage images directly from the terminal. The command-line interface provides powerful commands to build, run, and manage Docker containers. To install Docker on Linux, you can follow the instructions on the [Docker website](https://docs.docker.com/engine/install/).

### The Dockerfile

The Dockerfile is the blueprint for building the Docker image. It contains a series of commands that specify how to set up the environment, install necessary software, and configure the container. Below is an example of a Dockerfile that sets up an environment for scRNAseq analysis using R and Seurat.

Every docker image is built on top of a base image, this is the foundation on wich the other layer of the image are built. The base image can contain just the operating system or it can contain additional software and dependencies. In this case, we will use an Ubuntu base image, which is a popular choice for many Docker images due to its stability and wide range of available packages. 

You can specify which version of Ubuntu you want to use, in this case we will use the 20.04 version, which is a long-term support (LTS) release. This means that it will receive updates and support for an extended period, making it a reliable choice for building Docker images. If you do not specify a version, Docker will use the latest version available, which may not be compatible with all packages and dependencies.

```dockerfile
FROM ubuntu:20.04
# Set environment variables to avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive 
```

The base Ubuntu image is pretty bare, it does not contain any additional software or packages, just the operating system. To make it suitable for our tasks and for the installation of our required software, we need to install various packages. These packages include essential tools and libraries that are required for running R, Seurat, and Jupyter Lab. The following command installs a set of packages that are commonly used in data analysis and scientific computing.

```dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common dirmngr gpg curl build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libharfbuzz-dev \
    libfribidi-dev make cmake gfortran libxt-dev liblapack-dev libblas-dev \
    sudo wget zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev pandoc git && \
    rm -rf /var/lib/apt/lists/* # Clean up apt cache
```

Now we can install R, which is the programming language used for data analysis and visualization. The following commands add the CRAN GPG key and repository for R, update the package list, and install R. This ensures that we have the latest version of R installed in our Docker image.

```dockerfile
# Add the CRAN GPG key and repository for R
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    | tee /etc/apt/sources.list.d/cran-r.list
# Update again and install R
RUN apt update && apt install -y --no-install-recommends r-base
```

At the moment our docker image would have a functional R installation, but it could only be accesed via a terminal. To make it more user-friendly, we can install JupyterLab, which is a web-based interactive development environment for R and Python. JupyterLab allows you to create and share documents that contain live code, equations, visualizations, and narrative text. It is widely used in data science and scientific computing for its flexibility and ease of use.

To install JupyterLab, we first need to install Python and pip, the package manager for Python. We will also create a virtual environment to isolate the JupyterLab installation from the system Python environment. This is a good practice to avoid conflicts between different Python packages and versions.

With the last step we will install the IRkernel package, which allows JupyterLab to run R code. The `IRkernel::installspec(user = FALSE)` command registers the R kernel with JupyterLab, making it available for use in Jupyter notebooks.

```dockerfile
# Install JupyterLab
RUN apt update && apt install -y python3 python3-pip python3-venv
# create a virtual environment in which JupyterLab can be installed
RUN python3 -m venv /opt/venv
# Activate virtual environment and install JupyterLab
RUN /opt/venv/bin/pip install --upgrade pip && /opt/venv/bin/pip install jupyterlab
# Set the virtual environment as the default Python path
ENV PATH="/opt/venv/bin:$PATH"
# Make R visible to jupyter
RUN R -e "install.packages('IRkernel')" \
    R -e "IRkernel::installspec(user = FALSE)"
```

Now we can install the R packages that are required for the analysis. The following command installs the `Seurat` package, which is a popular R package for single-cell RNA sequencing data analysis. It also installs other useful packages such as `dplyr`, `ggplot2`, and `data.table` for data manipulation and visualization.

```dockerfile
# Install R packages
RUN R -e "install.packages(c('BiocManager', 'dplyr', 'ggplot2', 'data.table', 'future', 'cowplot', 'remotes', 'R.utils', 'rtracklayer', 'tinytex'))" 
RUN R -e "BiocManager::install(c('tidyverse', 'Seurat'))" 
```

A We want to set the user as a non-root user to avoid permission issues when running R and JupyterLab. This is a good practice to enhance security and avoid potential conflicts with file permissions. The following commands create a new user with the username `containeruser`, set the home directory, and switch to that user. (This part is not necessary and can be skipped, in case you want to run everything as root or in case it couses issues with your docker image)

```dockerfile
ARG USERNAME=containeruser

RUN useradd -u 1000 -m -s /bin/bash $USERNAME
USER $USERNAME
WORKDIR /home/$USERNAME
```

The last line in a Dockerfile is the command that will be executed when the container is started. 

In this case, we are setting the default shell to `/bin/bash` and starting JupyterLab with specific options. We are starting JupyterLab on port 8888, allowing access from any IP address without ne necessity of a token. This allows you to access JupyterLab from your web browser by going to `http://localhost:8888`. The `--allow-root` option allows JupyterLab to run as the root user, which is necessary when running inside a Docker container.

```dockerfile
# Set /bin/bash as the default shell
ENV SHELL=/bin/bash
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--ServerApp.allow_origin='*'", "--ServerApp.token=''"]
```

### Building the Docker image

To build the Docker image, you need to create the `Dockerfile` in your working directory. Then you can build the Docker image using the `docker build` command. This command takes the path to the directory containing the Dockerfile as an argument and builds the image according to the instructions in the Dockerfile.

The following command builds the Docker image and tags it with the name `scrnaseq_tutorial:latest`. The `-t` option specifies the name and tag of the image, and the `.` at the end indicates that the Dockerfile is in the current directory.

```sh
docker build -t scrnaseq_tutorial:latest .
```

It will take some time to build the image, as it needs to download and install all the packages and dependencies specified in the Dockerfile. Once the build is complete, you will see a message indicating that the image has been built successfully. Each layer of the image is cached, so if you rebuild the image without changing the Dockerfile, it will be much faster as it will reuse the cached layers.

### Running a Docker container

Now that the image has been built, you can run a Docker container from the it. A container is an instance of the Docker image that runs in an isolated environment. You can run multiple containers from the same image, each with its own isolated environment.

The command to run a Docker container is `docker run`. This command takes several options and arguments to specify how the container should be run. The most common options are:

| Option | Description |
|--------|-------------|
| `-it` | Run the container in interactive mode with a terminal |
| `--rm` | Automatically remove the container when it is stopped |
| `-p <host_port>:<container_port>` | Map a port on the host machine to a port in the container |
| `-v <host_dir>:<container_dir>` | Mount a directory from the host machine into the container, you can mount more that one directory by using multiple `-v` options |
| `--name <container_name>` | Assign a name to the container |
| `--network <network_name>` | Connect the container to a specific Docker network |
| `-d` | Run the container in detached mode, meaning it runs in the background and does not block the terminal |

In its simplest form, the command to run a Docker container is:

```sh
docker run -it ubuntu
```

This command runs a Docker container from the `ubuntu` image in interactive mode with a terminal. The `-it` option allows you to interact with the container's shell.

To run the container from the image we just built, you can use the following command:

```sh
docker run -it --rm -p 8888:8888 -v /path/of/working/directory:/sharedFolder scrnaseq_tutorial:latest
```
This command runs a Docker container from the `scrnaseq_tutorial:latest` image in interactive mode with a terminal. It maps port 8888 on the host machine to port 8888 in the container, allowing you to access JupyterLab from your web browser by goint to `http://localhost:8888`. It also mounts the `/path/of/working/directory` directory from the host machine into the `/sharedFolder` directory in the container, allowing you to access files from your working directory inside the container.

#### Run your analysis with one command
You can use docker containers to run your complete analysis. 

For organizational purposes you can mount two directories, one for the data and one for the results, this way you can keep them separate from the working directory. You can do this by using two `-v` options, one for each directory. 

With the `-d` option, you can run this container in detached mode, so that it runs in the background and does not block the terminal. 

By writing the name of a script after the name of the image you will be able to run it when the container is started. This is equivalent to the `CMD` command in the Dockerfile, but it will override it for this specific instance of the container.

(The `\` allows to split the command into multiple lines for better readability)

```sh
docker run -d --rm -p 8888:8888 \
    -v /path/of/dir/Data:/Data \
    -v /path/of/dir/Results:/Results \
    scrnaseq_tutorial:latest Script.R
```

### Saving changes to an image

If you download additional packages or make changes to the container, they will not be saved in the original image. To save the changes, you need to commit the container to a new image. This is done with the `docker commit` command. This command creates a new image from the changes made in the running container. It requires two arguments:

- `<container_id>`: The ID of the running container you want to commit.
- `<new_image_name>`: The name you want to give to the new image.

To see the ID of the running container, you can use the `docker ps` command, which lists all running containers along with their IDs.

```sh
docker ps
-> CONTAINER ID   IMAGE     COMMAND       CREATED          STATUS          PORTS     NAMES
-> a96481d90563   ubuntu    "/bin/bash"   55 seconds ago   Up 55 seconds             angry_shockley
docker commit a96481d90563 ubuntu:latest 
-> sha256:ed0d615f33cd2198cfd21d16ca8d79435b0e7fbe3c233c2aad07fb60cabd4b0b
```

You can commit the container to a new image with a specific name and tag. The tag is optional but recommended for versioning.

### Usefull Docker commands

This is a list of some useful Docker commands that can help you manage your Docker images and containers:

| Command | Description |
|---------|-------------|
| `docker images` | List all Docker images |
| `docker ps -a` | List all Docker containers, including stopped ones |
| `docker exec -it <container_id> /bin/bash` | Access a running Docker container's shell |
| `docker stop <container_id>` | Stop a running Docker container |
| `docker restart <container_id>` | Restart a stopped Docker container |
| `docker rm <container_id>` | Remove a specific Docker container |
| `docker rmi <image_id>` | Remove a specific Docker image |
| `docker container prune` | Remove all stopped Docker containers |
| `docker image prune` | Remove all unused Docker images |
| `docker volume prune` | Remove all unused Docker volumes |
| `docker network prune` | Remove all unused Docker networks |
| `docker system prune -a` | Remove all stopped containers, unused images, and unused networks |

## Preexisting Docker images

Many prebuilt Docker images are available on the Docker Hub repository, which can be used as a base for your own Dockerfile or run directly. This can save time and effort in setting up the environment, as these images often come with the necessary software and dependencies already installed.

If you try to run a container from an image that is not present on your system, Docker will automatically download it from the Docker Hub repository. You can also pull a specific image from the repository using the `docker pull` command.

The image build from the Dockerfile above is available on GitHub Container Registry, and can be pulled with the following command:

```sh
docker pull ghcr.io/maiolino-au/scrnaseq_tutorial:latest
```

Antoher one that might be usefull is the `satijalab/seurat` image, build by the Satija Lab, the developers of Seurat, and is available on [Docker Hub](https://hub.docker.com/r/satijalab/seurat). It contains the latest version of Seurat and all its dependencies. It runs on a terminal and does not contain JupyterLab, but it is possible to install it inside the container or to use this image as a base for your own Dockerfile. There is one available for each version of Seurat, so you can choose the one that best fits your needs. 

The following command pulls the Seurat 5.0.0 image from Docker Hub:

```sh
docker pull satijalab/seurat:5.0.0
```

## Credo: a better and faster way to build Docker images - WORK IN PROGRESS

[Credo](https://github.com/CREDOProject/core?tab=readme-ov-file) allows to build Docker images with specific versions of software, packages, and dependencies in a seemless and efficient way. 

We will use it for this project, but at the moment it is still in development and not all features are implemented. You do not have to build your Dockers with the exact versions specified in the papers, we will do it later with credo when it is ready.



# Singel Cell RNA Sequencing (scRNAseq) - WORK IN PROGRESS

what it is

## the data

Reads

allignemnt

## Dimensional reduction - WORK IN PROGRESS

### Linear methods - PCA - WORK IN PROGRESS

### Non-linear methods - UMAP, tSNE - WORK IN PROGRESS

## Clustering - WORK IN PROGRESS

## Cluster markers - WORK IN PROGRESS

## Pseudo Bulk - WORK IN PROGRESS

## Trajectory Inference - WORK IN PROGRESS

### Monocle - WORK IN PROGRESS

### RNA velocity - WORK IN PROGRESS




# Obtaining the data - WORK IN PROGRESS



## Downloading the data from GEO - WORK IN PROGRESS
The Gene Expression Omnibus (GEO) is a public database that stores high-throughput gene expression data, including single-cell RNA sequencing datasets. In this tutorial, we will use a dataset from GEO to demonstrate the analysis of scRNAseq data using the Seurat package in R.
The data can be downladed manually from the GEO site or with R.

### Manual download - WORK IN PROGRESS

### Download with R - WORK IN PROGRESS
The GEO database 
Expecially when working with multiple datasets it is faster to download the data directly from R. The dataset used in this tutorial is GSE150728, which contains scRNAseq data from human peripheral blood mononuclear cells (PBMCs) infected with SARS-CoV-2.

To dowload with R, you will need to install the `curl` and `R.utils` packages if they are not already installed.

```R
# Load required packages
if (!requireNamespace("curl", quietly = TRUE)) install.packages("curl")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")

library(curl)
library(R.utils)
```


```R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Provide a GEO id (ex. GSE150728)")
}
GEO_id <- args[1]

# URL and destination
ftp_tar <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    gsub(".{3}$", "", GEO_id), "nnn/",
    GEO_id, "/suppl/",
    GEO_id, "_RAW.tar"
)
dest_dir <- "/sharedFolder/Data"
tar_file <- file.path(dest_dir, "GSE150728_RAW.tar")

# Ensure directory exists
dir.create(dest_dir, showWarnings = FALSE)

# Download using curl with extended timeout
curl::curl_download(url = ftp_tar, destfile = tar_file, mode = "wb", handle = new_handle(timeout = 300))

# Extract the .tar archive
untar(tar_file, exdir = dest_dir)

# Unzip all .rds.gz files
gz_files <- list.files(dest_dir, pattern = "\\.gz$", full.names = TRUE)
for (f in gz_files) {
    message("Unzipping: ", f)
    gunzip(f, overwrite = TRUE, remove = TRUE)
}

# Optional cleanup
unlink(tar_file)
```



# Processing - WORK IN PROGRESS



## Load the data in R - WORK IN PROGRESS

The first step in the analysis of datas it is to load them in your working environment. But first we need to understant how they are structured.

After allignemnt, scRNA-seq data is typically stored in a sparse matrix. A full matrix is a table, with row and columns, in which every value is stored. This is not efficient for scRNA-seq data, as most of the value in the matrix will be zero (the majority of the genes are not expresset in a given cell). A sparse matrix is a data structure that only stores non zero values, greatly reducing the memory usage. A sparse matrix is divided in three files:

- `matrix.mtx.gz`: the sparse matrix in Matrix Market format, which contains the non-zero values of the matrix, along with their row and column indices.
- `features.tsv.gz`: a tab-separated file that contains the gene names and their corresponding indices in the sparse matrix. This file is used to map the gene names to their indices in the matrix.
- `barcodes.tsv.gz`: a tab-separated file that contains the cell barcodes and their corresponding indices in the sparse matrix. This file is used to map the cell barcodes to their indices in the matrix.

We therefore have two files with the names of rows and columns, and one with the value of the matrix and their coordinates. We can reassemble the sparse matrix in R.

### Load the spares matrix separately - WORK IN PROGRESS


```R
library(Matrix)
library(stats)

# Load matrix and metadata
mat <- readMM("/sharedFolder/Data/matrix.mtx.gz")
genes <- readLines("/sharedFolder/Data/features.tsv.gz")
barcodes <- readLines("/sharedFolder/Data/barcodes.tsv.gz")
```

### Load via Read10X - WORK IN PROGRESS

The `Read10X` function from the Seurat package can be used to read the sparse matrix and the metadata files in one go. This function is specifically designed to read 10X Genomics data, which is the format used for scRNA-seq data.

```R
raw_data <- Read10X(
    data.dir = "/sharedFolder/Data",
    gene.column = 1
)
```

Seurat obj

```R
sc_data <- CreateSeuratObject(
    counts = raw_data,
    min.cells = 3,
    min.features = 500,
    project = "scRNAseq_tutorial",
    names.delim = "-",
    names.field = 2
)
```

## Reads allignment - WORK IN PROGRESS - to move



```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

## Filter the cells - WORK IN PROGRESS



```R
# check whether you are getting the genes you want
grep("^mt-", rownames(sc_data@assays$RNA@data), value = T)
grep("Rps|Rpl|Mrpl|Mrps", rownames(sc_data@assays$RNA@data), value = T)
```


```R
# add the percentage of these mitochondrial or ribosomal genes to the meta.data
sc_data[["percent.mito"]] <- PercentageFeatureSet(
    object = sc_data,
    pattern = "^mt-"
)
sc_data[["percent.ribo"]] <- PercentageFeatureSet(
    object = sc_data,
    pattern = "Rps|Rpl|Mrpl|Mrps"
)
```


```R
# Subset the data based on the number of features and the percentage of mitochondrial and ribosomal genes
sc_data <- subset(
    x = sc_data,
    subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < 25 & percent.ribo < 40
)
```


## Normalization - WORK IN PROGRESS



```R
sc_data <- NormalizeData(
    sc_data,
    normalization.method = "LogNormalize",
    scale.factor = 1e6
)
```

## Variable features and Scaling - WORK IN PROGRESS



```R
# Find variable features
sc_data <- FindVariableFeatures(
    sc_data,
    selection.method = "mvp",
    nfeatures = 2000
)

# Scale the data
sc_data <- ScaleData(sc_data)
```

## Dimensionality Reduction - WORK IN PROGRESS



```R
# PCA
sc_data <- RunPCA(
    sc_data,
    npcs = 40
)
```


```R
# UMAP
sc_data_UMAP <- RunUMAP(sc_data, dims = 1:40)
```


```R
# Visualize UMAP
UMAP_plot <- DimPlot(sc_data_UMAP, reduction = "umap", label = TRUE, pt.size = 1) +
    ggtitle("UMAP of scRNAseq data")
print(plot)

# Save the UMAP plot
ggsave(
    filename = "UMAP_scRNAseq_data.png",
    plot = UMAP_plot,
    width = 1920, height = 1080, units = "px"
)
```

## Clustering - WORK IN PROGRESS

### Clusterization method 1 - WORK IN PROGRESS

```R
sc_data <- FindNeighbors(sc_data, dims = 1:40)
sc_data <- FindClusters(sc_data, resolution = 1)
```

### Clusterization method 2 - WORK IN PROGRESS

```R
sc_data <- FindNeighbors(dasc_datata, dims = 1:40)
sc_data <- FindClusters(sc_data, resolution = 1)
```

## Finding the most expressed markers - WORK IN PROGRESS



```R
markers <- FindAllMarkers(
    sc_data,
    only.pos = TRUE, # Only considers positive markers
    min.pct = 0.25, # Minimum percentage of cells expressing the gene
    logfc.threshold = 0.25, # Minimum log fold change
)
```

## Cell type annotation - WORK IN PROGRESS

### EnrichR - WORK IN PROGRESS

```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

### SingleR - WORK IN PROGRESS

```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

### GSEA - WORK IN PROGRESS

```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

### GSVA - WORK IN PROGRESS

```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

### GPTCellType - WORK IN PROGRESS

```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```



# Analysis - WORK IN PROGRESS



```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

## Differential Expression Analysis - WORK IN PROGRESS



```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```

## Visualization - WORK IN PROGRESS



```R
# mmmmmm
mmmm <- data.frame("m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
```
