## correlation matrix
## build gene matrix，use provided gene list
##step1: extract expression matrix of provided gene list(raw counts are fine) for each guide in the library, use median expression value
library(dplyr)
library(Seurat)
count <- as.matrix(GetAssayData(dat.var.seurat@assays$RNA, slot = "counts")) # 13254*3796 cells
pheno <- dat.var.seurat@meta.data[,c(5,6,11)]
idx <- row.names(count) %in%gene.list
count.gene.list <- count[idx,]
count.gene.list <- data.frame(pheno$feature_call,t(count.gene.list))
library(dplyr);library(stringr)
count.gene.median <- count.gene.list %>% group_by(pheno.feature_call) %>%
  summarise_all(median)
count.gene.median.pos <-
  count.gene.median %>%
  mutate(pos=as.numeric(word(pheno.feature_call,-2,sep=fixed('_')))) %>%
  arrange(pos)
genes_to_keep <- colSums(count.gene.median.pos)!=0
count.gene.median.pos.rm0 <- count.gene.median.pos[,genes_to_keep]
## step2: calculate the cor coefficient for each guide-pair
count.gene.median.pos.rm0.tr <- t(count.gene.median.pos.rm0)
dat.cormat <- round(cor(count.gene.median.pos.rm0.tr),2)
View(head(dat.cormat))
## step 3: draw heatmap
library(reshape2)
melted_dat_cormat <- melt(dat.cormat)
#heatmap all
ggplot(data = melted_dat_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()
# Heatmap upper
ggplot(data = melted_dat_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))+
  coord_fixed()
