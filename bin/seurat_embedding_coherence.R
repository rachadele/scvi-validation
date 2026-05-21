#!/usr/bin/env Rscript
#
# seurat_embedding_coherence.R
#
# Computes pairwise cosine similarity between ref and query cell-type centroids
# in Seurat PCA space (query projected into ref PCA via pcaproject).
# Mirrors embedding_coherence.py but uses the Seurat integration embedding
# instead of scVI latent space.
#
# Outputs (per ref_key):
#   {query_name}_{ref_name}.{key}.seurat_coherence.tsv  -- matrix (ref types x query types)
#   {query_name}_{ref_name}.{key}.seurat_coherence.png  -- heatmap

library(Seurat)
library(argparse)

cosine_similarity_matrix <- function(A, B) {
    # A: n x d, B: m x d  ->  returns n x m cosine similarity matrix
    A_norm <- A / sqrt(rowSums(A^2))
    B_norm <- B / sqrt(rowSums(B^2))
    return(A_norm %*% t(B_norm))
}

compute_centroids <- function(embeddings, labels) {
    unique_labels <- sort(unique(labels))
    centroids     <- do.call(rbind, lapply(unique_labels, function(lbl) {
        colMeans(embeddings[labels == lbl, , drop = FALSE])
    }))
    rownames(centroids) <- unique_labels
    return(centroids)
}

parser <- argparse::ArgumentParser(description = "Seurat embedding centroid cosine similarity.")
parser$add_argument("--ref_path",   type = "character", required = TRUE)
parser$add_argument("--query_path", type = "character", required = TRUE)
parser$add_argument("--ref_keys",   type = "character", nargs = "+",
                    default = c("subclass", "class", "family", "global"))
parser$add_argument("--dims",       type = "integer",   default = 50)
parser$add_argument("--seed",       type = "integer",   default = 42)

args       <- parser$parse_args()
set.seed(args$seed)

ref_name   <- gsub(".rds$", "", basename(args$ref_path))
query_name <- gsub(".rds$", "", basename(args$query_path))
prefix     <- paste0(query_name, "_", ref_name)

message("Loading ref: ", args$ref_path)
ref   <- readRDS(args$ref_path)
message("Loading query: ", args$query_path)
query <- readRDS(args$query_path)

# Project query into ref PCA space
message("Finding transfer anchors (pcaproject)...")
anchors <- FindTransferAnchors(
    reference           = ref,
    query               = query,
    normalization.method = ifelse("SCT" %in% names(ref@assays), "SCT", "LogNormalize"),
    reference.reduction = "pca",
    reduction           = "pcaproject",
    dims                = 1:args$dims
)

# Use IntegrateEmbeddings to project query into ref PCA space.
# This is the embedding-projection half of MapQuery, without TransferData
# (label transfer), which is what triggers the k.weight error.
message("Integrating embeddings (projecting query into ref PCA space)...")
query <- tryCatch({
    IntegrateEmbeddings(
        anchorset          = anchors,
        reference          = ref,
        query              = query,
        new.reduction.name = "ref.pca"
    )
}, error = function(e) {
    message("IntegrateEmbeddings error: ", e$message)
    if (grepl("less than", e$message)) {
        match    <- regmatches(e$message, regexec("less than ([0-9]+)", e$message))
        k.weight <- as.integer(match[[1]][2]) - 1L
        message("Retrying with k.weight = ", k.weight)
        IntegrateEmbeddings(
            anchorset          = anchors,
            reference          = ref,
            query              = query,
            new.reduction.name = "ref.pca",
            k.weight           = k.weight
        )
    } else {
        stop(e)
    }
})

ref_embeddings   <- Embeddings(ref,   "pca")
query_embeddings <- Embeddings(query, "ref.pca")

for (key in args$ref_keys) {
    ref_col   <- ref@meta.data[[key]]
    query_col <- query@meta.data[[key]]

    if (is.null(ref_col) || is.null(query_col)) {
        message("Skipping ", key, ": not in both ref and query metadata")
        next
    }

    ref_labels   <- as.character(ref_col)
    query_labels <- as.character(query_col)

    ref_centroids   <- compute_centroids(ref_embeddings,   ref_labels)
    query_centroids <- compute_centroids(query_embeddings, query_labels)

    sim_matrix <- cosine_similarity_matrix(ref_centroids, query_centroids)

    tsv_file <- paste0(prefix, ".", key, ".seurat_coherence.tsv")
    df       <- as.data.frame(sim_matrix)
    write.table(df, file = tsv_file, sep = "\t", quote = FALSE, col.names = NA)
    message("Wrote ", tsv_file)
}

message("Done.")
