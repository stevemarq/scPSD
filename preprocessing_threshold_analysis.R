suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(Seurat)
})


matrix <- "/Users/stevem/Desktop/bioinform/scPSD/datasets/matrix.mtx.gz"
features <- "/Users/stevem/Desktop/bioinform/scPSD/datasets/genes.tsv.gz"
cells <- "/Users/stevem/Desktop/bioinform/scPSD/datasets/barcodes.tsv.gz"

data <- ReadMtx(mtx = matrix, features = features, cells = cells)
seurat_data <- CreateSeuratObject(counts = data)
seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
VlnPlot(seurat_data, features <- c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
