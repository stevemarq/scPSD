library(Seurat)
library(dplyr)
#BiocManager::install("Linnorm") <-- so somehting with this somehow
library(Linnorm)
library(edgeR)
library(scone)

## add some code here so that the user inserts the wdir
data <- ReadMtx(
  mtx = "Desktop/bioinform/scPSD/GSE157829_RAW/GSM4775588_C1matrix.mtx.gz", 
  features = "Desktop/bioinform/scPSD/GSE157829_RAW/GSM4775588_C1genes.tsv.gz", 
  cells = "Desktop/bioinform/scPSD/GSE157829_RAW/GSM4775588_C1barcodes.tsv.gz"
)


seurat_data <- CreateSeuratObject(counts = data)
seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
target <- subset(seurat_data, subset = nCount_RNA > 200 & percent.mt < 30 
                 & nFeature_RNA > 0)
barcodes <- rownames(target)
features <- colnames(target)

data <- GetAssayData(target, layer = "counts")
data <- as.matrix(data)

processed <- Linnorm(data)







