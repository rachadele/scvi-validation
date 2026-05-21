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

parser <- argparse::ArgumentParser(description = "Convert H5AD to RDS (Seurat) for reference.")
parser$add_argument("--h5ad_file",            type = "character", required = TRUE)
parser$add_argument("--normalization_method", type = "character", default = "SCT")
parser$add_argument("--dims",                 type = "integer",   default = 50)
parser$add_argument("--batch_key",            type = "character", default = "dataset_title")
parser$add_argument("--nfeatures",            type = "integer",   default = 2000)
parser$add_argument("--batch_correct",        action = "store_true", default = TRUE)
parser$add_argument("--k.anchor",             type = "integer",   default = 5)
parser$add_argument("--k.weight",             type = "integer",   default = 30)
parser$add_argument("--k.score",              type = "integer",   default = 30)

args                 <- parser$parse_args()
normalization_method <- args$normalization_method
dims                 <- args$dims
batch_key            <- args$batch_key
n_features           <- args$nfeatures
h5ad_file            <- args$h5ad_file
batch_correct        <- args$batch_correct
k.anchor             <- args$k.anchor
k.weight             <- args$k.weight
k.score              <- args$k.score

sceasy_seurat <- sceasy::convertFormat(h5ad_file, from = "anndata", to = "seurat")

if ("feature_id" %in% colnames(sceasy_seurat@assays$RNA[[]])) {
    sceasy_seurat <- rename_features(sceasy_seurat, column_name = "feature_id")
}

if (batch_correct) {
    sceasy_list <- SplitObject(sceasy_seurat, split.by = batch_key)

    if (length(sceasy_list) == 1) {
        message("Only one batch found. Skipping integration.")
        batch_correct <- FALSE
    } else {
        sceasy_list <- Filter(function(x) ncol(x) > dims, sceasy_list)

        if (normalization_method == "SCT") {
            assay <- "SCT"
            sceasy_list <- lapply(sceasy_list, function(x) SCTransform(x, verbose = FALSE, variable.features.n = n_features))
            features    <- SelectIntegrationFeatures(sceasy_list)
            sceasy_list <- PrepSCTIntegration(sceasy_list, anchor.features = features)
        } else if (normalization_method == "LogNormalize") {
            assay <- "RNA"
            sceasy_list <- lapply(sceasy_list, function(x) {
                x <- NormalizeData(x)
                x <- FindVariableFeatures(x, nfeatures = n_features)
                x <- ScaleData(x)
                return(x)
            })
            features <- SelectIntegrationFeatures(sceasy_list)
        }

        anchors <- FindIntegrationAnchors(
            object.list          = sceasy_list,
            dims                 = 1:dims,
            k.anchor             = k.anchor,
            k.score              = k.score,
            normalization.method = normalization_method,
            anchor.features      = features
        )

        anchors_per_pair <- table(anchors@anchors$dataset1, anchors@anchors$dataset2)
        min_anchors      <- min(anchors_per_pair[anchors_per_pair > 0])
        k.anchor         <- min(min_anchors, k.anchor)

        sceasy_seurat <- IntegrateData(
            anchorset            = anchors,
            dims                 = 1:dims,
            k.weight             = k.weight,
            normalization.method = normalization_method,
            new.assay.name       = "integrated"
        )
        assay <- "integrated"
        if (normalization_method == "LogNormalize") {
            sceasy_seurat <- sceasy_seurat %>% ScaleData(assay = assay)
        }
    }
}

if (!batch_correct) {
    if (normalization_method == "SCT") {
        sceasy_seurat <- SCTransform(sceasy_seurat, verbose = FALSE, variable.features.n = n_features)
        assay <- "SCT"
    } else if (normalization_method == "LogNormalize") {
        sceasy_seurat <- NormalizeData(sceasy_seurat, normalization.method = normalization_method)
        sceasy_seurat <- FindVariableFeatures(sceasy_seurat, nfeatures = n_features)
        sceasy_seurat <- ScaleData(sceasy_seurat)
        assay <- "RNA"
    }
}

sceasy_seurat <- RunPCA(sceasy_seurat, npcs = dims, assay = assay)

out_file <- gsub(".h5ad", ".rds", basename(h5ad_file))
saveRDS(sceasy_seurat, file = out_file)
message(paste("Saved to", out_file))
