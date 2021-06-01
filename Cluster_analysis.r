


### add guide info, known functional region info to individual cells
library(dplyr)
library(stringr)
crispr_gene_feature_ref <- read.csv("YOUR_REF_FILE_PATH/crispr_gene_feature_ref.csv", stringsAsFactors=FALSE)
temp<-crispr_gene_feature_ref %>% mutate(pos=as.numeric(word(id,-2,sep=fixed('_'))))%>%select(id,pos)
## region known pos example
temp <- within(temp, category[pos <= 378] <- "KMT_Nterm")
## create protospacer info and meta info, note final meta info rownames must be cellbarcodes identifier
protospacer_calls_per_cell_dmso <- read.csv("/Volumes/CWChen_Grp/Novogene/X202SC20021687-Z01-F001/raw_data/DMSO_Dot1l_0331/outs/crispr_analysis/protospacer_calls_per_cell.csv", stringsAsFactors=FALSE)

save(temp,protospacer_calls_per_cell,file="annofiles.RData")


### guide clustering analysis based on transcriptome info
library(Seurat);library(dplyr)
dat.seurat <- readRDS(file = "./dat.rd") # corrected dmso
meta.info <- dat.seurat@meta.data
load("./annofiles.RData")
meta.info.new <- merge(meta.info,protospacer_calls_per_cell,all.x = T,by.x = "row.names",by.y="cell_barcode")
meta.info.new$cell_barcode <- meta.info.new$Row.names
row.names(meta.info.new) <- meta.info.new$Row.names  ## 14252*14,  NAs
meta.info.new <- meta.info.new[,-1]
## identify cells with zero guide and more than one guide
meta.info.new$num_features <- ifelse(is.na(meta.info.new$num_features), 0, meta.info.new$num_features)
meta.info.new$feature_call <- ifelse(is.na(meta.info.new$feature_call), "Empty", meta.info.new$feature_call)
meta.info.new$num_umis <- ifelse(is.na(meta.info.new$num_umis), 0, meta.info.new$num_umis)
idx <- grep(pattern = "\\|",x = meta.info.new$feature_call)
meta.info.new$feature_call[idx] <- "more_than_one_guide"
meta.info.new$category <- sapply(meta.info.new$feature_call, function(x)add_category(x,temp))
View(as.data.frame(table(meta.info.new$category)))
dat.seurat@meta.data <- meta.info.new
## subset by excluding empty and more than one guide cells and pos_ctrl
meta.info.sub <- meta.info.new %>%
  select(orig.ident,nCount_RNA,nFeature_RNA,percent.mt,feature_call,category,num_umis)
dat.seurat@meta.data <- meta.info.sub    ## 6704*7
Idents(dat.seurat) <- "category"
dat.seurat.sub <- subset(dat.seurat, idents = c("AF9","Neg_Ctrl" ,"KMT_SAM","KMT_Nterm", "R1", "R2","UKND"))   ## replace the catgory name with your own ones

### run pca using provided gene list gene.list
dat.var.seurat <- RunPCA(dat.seurat.sub, features = gene.list)
## umap plot
dat.var.seurat <- RunUMAP(dat.var.seurat, dims = 1:30,umap.method="umap-learn")
dat.var.seurat <- FindNeighbors(dat.var.seurat, dims = 1:30)
dat.var.seurat <- FindClusters(dat.var.seurat, resolution = 0.5)
pdf("UMAP_plot.pdf",width = 9, height = 9)
#DimPlot(dat.var.seurat, reduction = "umap", label = T)
DimPlot(dat.var.seurat,reduction = "umap", split.by = "category", group.by = "seurat_clusters",pt.size=0.5, label = T, label.size=6)
dev.off()




## heatmap correaltion analysis
