#BiocManager::install("scran") # <-- so something with this somehow
.rs.restartR()
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(scone)
  library(Linnorm)
  library(scran)
  library(Seurat)
  library(argparse) 
  })

print("Necessary libraries loaded")


load_data <- function(matrix, features, cells){
  
  data <- ReadMtx(mtx = matrix, features = features, cells = cells)
  return (data)
}


preprocess <- function(data, 
                       norm = c("tmm", "cpm", "scone", "linnorm", "scran", 'surat'),
                       nCRNA_threshold, nFRNA_threshold, pmt){
  
  seurat_data <- CreateSeuratObject(counts = data)
  seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
  target <- subset(seurat_data,
                   subset = nCount_RNA > nCRNA_threshold 
                   & percent.mt < pmt 
                   & nFeature_RNA > nFRNA_threshold)
  barcodes <- rownames(target)
  features <- colnames(target)
  
  data <- GetAssayData(target, layer = "counts")
  data <- as.matrix(data)
  
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
    errorCondition("NotImplemented")
  }else if (norm == 'linnorm'){
    norm_data <- Linnorm.Norm(data, output = "Raw")
  }else if (norm == 'scran'){
    
  }
  
  
  metrics <- list(norm_data, barcodes, features)
  return (metrics)
}




parser <- ArgumentParser(description = "Loading in GEO dataset and 
                         preprocessing before clustering.")

# dataset loader
parser$add_argument("-m", "--mtx", type = "character", help = " The genetic .mtx.gz file" )
parser$add_argument("-b", "--barcode", type = "character", 
                    help = " The cell barcodes in a .tsv.gz file" )
parser$add_argument("-g", "--genes", type = "character",
                    help = "The genes in a .tsv.gz file" )


# data clean up 
parser$add_argument("-cr", "--cRNA", type = "integer", default = 200, 
                    help = "Min threshold for the number of RNA counts per cells")
parser$add_argument("-fr", "--fRNA", type = "integer", default = 0, 
                    help = "Min threshold for the number of genes expressed per cells")
parser$add_argument("-p", '--pmt', type = "integer", default = 30,
                    help = "Max threshold for percent of mitochondrial content per cell")
parser$add_argument("-n", type = "character", default = "Raw",
                    help = "Type of matrix normalization")



args <- parser$parse_args()


# location of datasets
script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", commandArgs[grep("--file=", script_path)])
script_dir <- dirname(normalizePath(script_path))

# "/Users/stevem/Desktop/bioinform/scPSD"
mtx <- paste(script_dir, args$mtx, sep="/")
barcode <- paste(script_dir, args$barcode, sep="/")
features <- paste(script_dir, args$genes, sep="/")
nCRNA_threshold <- args$cRNA
nFRNA_threshold <- args$fRNA
pmt <- args$pmt
norm <- args$norm

data_matrix <- load_data(mtx, features, barcodes)
new_data_matrix <- preprocess(data_matrix, norm, nCRNA_threshold, nFRNA_threshold, pmt)

data("Islam2011")
dim(data)
y_norm<- Linnorm.Norm(data, output='Raw')
