#BiocManager::install("scran") # <-- so something with this somehow

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(edgeR)
  library(scone)
  library(Linnorm)
  library(scran)
  library(Seurat)
  library(argparse)
  library(celldex)
  library(SingleR)
  })
#BiocManager::install(c("SingleR", "celldex"))

print("Necessary libraries loaded")



get_cell_types <- function(data){
    ref <- celldex::BlueprintEncodeData() # uses this one because we know we are looking at immune cells
    pred <- SingleR(test=data, ref=ref, labels=ref$label.main)
    classes <- pred$labels
    classes_df <- as.data.frame(classes)
    return (classes_df)
}


load_data <- function(matrix, features, cells){
  
  for (f in c(matrix, features, cells)){
    if (!file.exists(f)){
      stop(paste(f, "file does not exist", sep=" "))
      }
  }
  
  data <- ReadMtx(mtx = matrix, features = features, cells = cells)
  return (data)
}


preprocess <- function(data, 
                       norm = c("raw", "tmm", "cpm", "scone", "linnorm", "scran", 'surat'),
                       nCRNA_threshold, nFRNA_threshold, pmt){
  
  
  
  seurat_data <- CreateSeuratObject(counts = data)
  seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
  target <- subset(seurat_data,
                   subset = nCount_RNA < nCRNA_threshold 
                   & percent.mt < pmt 
                   & nFeature_RNA < nFRNA_threshold)
  
  data <- GetAssayData(target, layer = "counts")
  features <- rownames(data)
  cell_types <- get_cell_types(data)
  
  # check normalization types
  norm <- tolower(norm)    
  norm <- match.arg(norm)

  if (norm == 'tmm'){
    dge <- DGEList(counts = data)
    dge <- calcNormFactors(dge, method = 'TMM')
    norm_data <- edgeR::cpm(dge, log=FALSE)
    
  }else if ( norm == 'cpm'){
    norm_data <- edgeR::cpm(data)
  }else if (norm == 'scone'){
    stop("NotImplemented")
  }else if (norm == 'linnorm'){
    norm_data <- Linnorm.Norm(data, output = "Raw")
  }else if (norm == 'scran'){
    clusters <- quickCluster(data)
    sce <- computeSumFactors(data, clusters=clusters)
    norm_data <- logNormCounts(sce)
  } else {
    norm_data <- data
  }



  metrics <- list(mtx = norm_data,
                  cells = cell_types,
                  features = features)
  return (metrics)
}




parser <- ArgumentParser(description = "Loading in GEO dataset and 
                         preprocessing before clustering.")

# dataset loader
parser$add_argument("-m", "--mtx", type = "character", help = " The genetic file '[dataset_name]/[filename].mtx.gz'" )
parser$add_argument("-b", "--barcode", type = "character", 
                    help = " The cell barcodes in a '[dataset_name]/[filename].tsv.gz' form" )
parser$add_argument("-g", "--genes", type = "character",
                    help = "The genes in a '[dataset_name]/[filename].tsv.gz' form" )


# data clean up 

parser$add_argument("-cr", "--cRNA", type = "integer", default = 2000, 
                    help = "Max threshold for the number of RNA counts per cells")
parser$add_argument("-fr", "--fRNA", type = "integer", default = 2000, 
                    help = "Max threshold for the number of genes expressed per cells")
parser$add_argument("-p", '--pmt', type = "integer", default = 30,
                    help = "Max threshold for percent of mitochondrial content per cell")
parser$add_argument("-n", "--norm", type = "character", default = "Raw",
                    help = "Type of matrix normalization")



args <- parser$parse_args()

# location of datasets
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_path[grep("--file=", script_path)])
script_dir <- dirname(normalizePath(script_path))

# script_dir == "/Users/stevem/Desktop/bioinform/scPSD"
mtx <- paste(script_dir, args$mtx, sep="/")
barcodes <- paste(script_dir, args$barcode, sep="/")
features <- paste(script_dir, args$genes, sep="/")
nCRNA_threshold <- args$cRNA
nFRNA_threshold <- args$fRNA
pmt <- args$pmt
norm <- args$norm

data_matrix <- load_data(mtx, features, barcodes)
print("Dataset Loading Complete")
new_data_matrix <- preprocess(data_matrix, norm, nCRNA_threshold, nFRNA_threshold, pmt)
print("Data Normalization Complete")

processed_data_file <- paste(script_dir, "processed_data", sep="/")
prep_out <- paste(processed_data_file, "prep_out.mtx", sep="/")
cells_out <- paste(processed_data_file, "cells_out.csv", sep="/")
feats_out <- paste(processed_data_file, "feats_out.csv", sep="/")


ifelse(!dir.exists(processed_data_file),
       dir.create(processed_data_file),
       "Folder exists already")


new_sparse_mtx <- as(new_data_matrix$mtx, "dgCMatrix")
writeMM(new_sparse_mtx, file = prep_out)
write.csv(new_data_matrix$cells, file = cells_out, row.names = FALSE)
write.csv(new_data_matrix$features, file = feats_out, row.names = FALSE)


