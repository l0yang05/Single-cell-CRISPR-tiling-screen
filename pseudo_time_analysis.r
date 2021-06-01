library(monocle)## version 2.14.0
library(Seurat)
               setwd("/Users/luyang/Documents/Projects/crispr_10x/lib610_corrected_DMSO")
# dat.var.seurat from umap clustering analysis
count <- as.matrix(GetAssayData(dat.var.seurat@assays$RNA, slot = "counts")) # 13254*3796 cells
## create featureData object
gene_anno <- data.frame(gene_short_name=dat.var.seurat@assays$RNA@counts@Dimnames[[1]],anno=array(dim=N))
row.names(gene_anno) <- gene_anno[,1]
fd <- new("AnnotatedDataFrame", data = gene_anno)
## create phenoData object
pheno <- dat.var.seurat@meta.data[,c(5,6,11)] # select anno info columns based on needs
pd <- new("AnnotatedDataFrame", data = pheno)
dat.mono <- newCellDataSet(count,phenoData = pd, featureData = fd, expressionFamily=negbinomial.size())
dat.mono <- estimateSizeFactors(dat.mono)
dat.mono <- estimateDispersions(dat.mono)

### Trajectory step 1: choose genes that define a cell's progress,
## odering genes selection (use seurat markers forthis test run )
library(dplyr)
dat.mono <- setOrderingFilter(dat.mono, ordering_genes = gene.list)
## Trajectory step 2: reduce data dimensionality
dat.mono <- reduceDimension(dat.mono, max_components = 2, method = 'DDRTree')

## Trajectory step 3: order cells along the trajectory
dat.mono <- orderCells(dat.mono)
plot_cell_trajectory(dat.mono, color_by = "seurat_clusters",show_branch_points=F)+
  facet_wrap(~seurat_clusters, nrow = 1)
plot_cell_trajectory(dat.mono, color_by = "category",show_branch_points=F,cell_size = 0.8)+
  facet_wrap(~category, nrow = 1)

plot_cell_trajectory(dat.mono, color_by = "Pseudotime",show_branch_points=F,cell_size = 0.8)+
  facet_wrap(~category, nrow = 1)

plot_genes_jitter(dat.mono["Meis1",], grouping = "category", min_expr = 0.1)
plot_genes_jitter(dat.mono["Hoxa9",], grouping = "category", min_expr = 0.1)
length(unique(pData(dat.mono)$State)) ## linear, only one state

## export monocle object
mono.anno <- as(dat.mono@phenoData, "data.frame")
mono.count <- exprs(dat.mono)
## create Seurat object
dat.mono.seurat<- CreateSeuratObject(counts = mono.count, project = "from_monocle",meta.data=mono.anno)
