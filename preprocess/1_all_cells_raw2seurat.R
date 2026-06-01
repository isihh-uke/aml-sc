#Sketch integration rerun 
#Seurat V5 
#zheng
#"Thu Dec 14 13:50:16 2023"
#15 sub-libs 777,530 cells

# packages and functions ---------------------------------------------------------------

library(Seurat) 
library(SeuratObject)
library(BPCells)

library(dplyr)
library(purrr)
library(ggrepel)
library(patchwork)
library(stringr) 
library(magrittr)

library(cowplot) 
library(ggrastr) 
library(ggplot2)

library(future)

# set this option when analyzing large datasets
options(future.globals.maxSize = 256*1024^3)
options(Seurat.object.assay.version = "v5")

future::plan(strategy = "multicore", workers = 32)

setwd("/path/to/your/wd")
wd_path <- getwd()
mat_path <- "/path/to/your/matrix"

# read matrix data --------------------------------------------------------
all_mat_h5ad <- BPCells::open_matrix_anndata_hdf5(paste0(mat_path,'/all_cells_raw.h5ad')
 
# Write the matrix to a directory
BPCells::write_matrix_dir(mat = all_mat_h5ad,
                          dir = paste0(wd_path, '/bpcell_mat'))

#read mat which is the output of BPCell
all_mat <- BPCells::open_matrix_dir(dir = paste0(wd_path, '/bpcell_mat'))

# read in cell metadata
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

# create seurat object
all_cells <- Seurat::CreateSeuratObject(counts = all_mat, meta.data = cell_meta)

rm(all_mat_h5ad)
rm(all_mat)
rm(cell_meta)

# QC ----------------------------------------------------------------------
# Befroe QC: 778438  cells
all_cells[["percent_mt"]] <- PercentageFeatureSet(all_cells, pattern = "^MT-")

colnames(all_cells@meta.data)[c(7,8)] <- c("nFeature_RNA", "nCount_RNA")

# Feature_rast are in-house developed function. 
# It is similar to Seurat::DimPlot

VlnPlot(all_cells, pt.size = 0.02, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, raster = T)

Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "nFeature_RNA", g = "percent_mt", color_grd = "grd", noaxis = F, axis.number = T)

Feature_rast(all_cells, d1 = "nCount_RNA", d2 = "percent_mt", g = "nFeature_RNA", color_grd = "grd", noaxis = F, axis.number = T)

# Subset Seurat
all_cells %<>% subset(nFeature_RNA < 10000 & nCount_RNA < 100000 & percent_mt < 50)
# After first filter: 778371 cells

# Subset Seurat
all_cells %<>% subset(nFeature_RNA < 8500 & nCount_RNA < 60000 & percent_mt < 25)
# After second filter: 778071 cells 

# Normalize
all_cells <- NormalizeData(all_cells)

# Split assay into layers
all_cells[["RNA"]] <- split(all_cells[["RNA"]], f = all_cells$donor)

# switch to analyzing the full dataset (on-disk)
DefaultAssay(all_cells) <- "RNA"

# Sketch representative cells ------------------------------------------------------
# Sketch 10k cells from total dataset
all_cells <- SketchData(all_cells, assay = "RNA", ncells = 10000, 
                        sketched.assay = "sketch", method = "LeverageScore")

DefaultAssay(all_cells) <- "sketch"

# Find HVGs                                                
all_cells %<>%  FindVariableFeatures(verbose = T, 
                                     selection.method = 'vst', 
                                     nfeatures = 2000)

top10 <- head(VariableFeatures(all_cells), 10)
plot1 <- VariableFeaturePlot(all_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# Dimentional reduciton ---------------------------------------------------------------
all_cells %<>% ScaleData(verbose = T) %>% RunPCA(verbose = T)

VizDimLoadings(all_cells, dims = 1:30, reduction = "pca")
Seurat::ElbowPlot(all_cells, ndims = 30, reduction = 'pca')
all_cells %<>% JackStraw(num.replicate = 100, dims = 30) %>% ScoreJackStraw(dims = 1:30)
JackStrawPlot(all_cells, dims = 1:30)

# Harmony integration
all_cells <- IntegrateLayers(object = all_cells, method = HarmonyIntegration,
                             orig.reduction = "pca", new.reduction = "harmony",
                             verbose = TRUE, assay = 'sketch')
# Cluster the integrated data
all_cells <- FindNeighbors(all_cells, reduction = "harmony", dims = 1:14)
all_cells <- FindClusters(all_cells, resolution = 2, cluster.name = "harmony_clusters")
all_cells <- RunUMAP(all_cells, reduction = "harmony", dims = 1:14, verbose = T)

# Project-integrate the full datasets 
all_cells <- ProjectIntegration(object = all_cells, sketched.assay = "sketch", 
                                assay = "RNA", reduction = "harmony")

all_cells <- ProjectData(object = all_cells, sketched.assay = "sketch", assay = "RNA",
                         sketched.reduction = "harmony.full", 
                         full.reduction = "harmony.full", 
                         dims = 1:14, 
                         refdata = list(harmony_clusters_full = "harmony_clusters")
                        )

all_cells <- RunUMAP(all_cells, reduction = "harmony.full", 
                     dims = 1:14, reduction.name = "umap.full", 
                     reduction.key = "UMAPfull_")

saveRDS(object = all_cells, file = paste0(wd_path, "/aml_all_cells_harmony.Rds"))