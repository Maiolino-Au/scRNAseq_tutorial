# Load required packages
if (!requireNamespace("curl", quietly = TRUE)) install.packages("curl")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
library(curl)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Provide a GEO id (ex. GSE150728)")
}
GEO_id <- args[1]

# URL and destination
ftp_tar <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", gsub(".{3}$", "", GEO_id), "nnn/", GEO_id, "/suppl/", GEO_id, "_RAW.tar")
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
