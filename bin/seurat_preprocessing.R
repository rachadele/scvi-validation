library(Seurat)
library(reticulate)
use_condaenv("/home/rschwartz/anaconda3/envs/r4.3/")
library(sceasy)
library(argparse)
library(tidyr)
options(future.globals.maxSize = 5 * 1024^3)  # 5 GB

rename_features <- function(seurat_obj, column_name) {
    counts <- seurat_obj[["RNA"]]@counts
    data   <- seurat_obj[["RNA"]]@data

    feature_meta <- seurat_obj[["RNA"]][[]]
    feature_meta[["orig_features"]] <- rownames(feature_meta)
    rownames(feature_meta) <- feature_meta[[column_name]]

    rownames(counts) <- rownames(feature_meta)
    rownames(data)   <- rownames(feature_meta)

    newRNA <- CreateAssayObject(counts = counts)
    newRNA$data <- data
    newRNA[[]]  <- feature_meta

    seurat_obj[["RNA"]] <- newRNA
    DefaultAssay(seurat_obj) <- "RNA"
    return(seurat_obj)
}

parser <- argparse::ArgumentParser(description = "Convert H5AD to RDS (Seurat) for query.")
parser$add_argument("--h5ad_file",            type = "character", required = TRUE)
parser$add_argument("--normalization_method", type = "character", default = "SCT")
parser$add_argument("--dims",                 type = "integer",   default = 50)
parser$add_argument("--nfeatures",            type = "integer",   default = 2000)

args               <- parser$parse_args()
h5ad_file          <- args$h5ad_file
normalization_method <- args$normalization_method
dims               <- args$dims
nfeatures          <- args$nfeatures

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from = "anndata", to = "seurat")

if ("feature_id" %in% colnames(sceasy_seurat@assays$RNA[[]])) {
    sceasy_seurat <- rename_features(sceasy_seurat, column_name = "feature_id")
}

if (normalization_method == "LogNormalize") {
    sceasy_seurat <- sceasy_seurat %>%
        NormalizeData(normalization.method = normalization_method) %>%
        FindVariableFeatures(nfeatures = nfeatures) %>%
        ScaleData() %>%
        RunPCA(npcs = dims)
} else if (normalization_method == "SCT") {
    sceasy_seurat <- sceasy_seurat %>%
        SCTransform(verbose = TRUE, variable.features.n = nfeatures) %>%
        RunPCA(npcs = dims, assay = "SCT")
} else {
    stop("Normalization method not recognized.")
}

out_file <- gsub(".h5ad", ".rds", basename(h5ad_file))
saveRDS(sceasy_seurat, file = out_file)
message(paste("Saved to", out_file))
