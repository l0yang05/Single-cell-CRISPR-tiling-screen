library(stringr)
# extract guide position, code may be modified based on your own datatable format
pos <- as.numeric(word(colnames(dat.cormat),-2,sep=fixed('_')))
dat.sm <- sapply(1:nrow(dat.cormat), function(i)smoothen.cor(rowdat=dat.cormat[i,], pos=pos))
colnames(test) <- pos
dat.sm.1 <- sapply(nrow(dat.sm), function(i)smoothen.cor(rowdat=dat.sm[i,], pos=pos))
colnames(dat.sm.1) <- row.names(dat.sm.1)
saveRDS(dat.sm.1,file="datCormat_smoothenByrow_cor.rd")
## example smoothen plot with DOT1L data
## upper matrix only
upper_tri_dmso_top100_row_cor <- get_upper_tri(dat.sm.1)
melted_dmso_top100_row_upper_cormat <- melt(upper_tri_dmso_top100_row_cor, na.rm = TRUE)
ggplot(data = melted_dmso_top100_row_upper_cormat, aes(Var2, Var1, fill = value))+geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 1,
                                   size = 8, hjust = 1))+
  coord_fixed()+
  geom_vline(xintercept = 362,color="blue", linetype="dotted")+
  geom_vline(xintercept = 938,color="blue", linetype="dotted")+
  geom_vline(xintercept = 1375,color="red", linetype="dotted") +
  geom_vline(xintercept = 1643,color="red", linetype="dotted") +
  geom_vline(xintercept = 1679,color="green", linetype="dotted") +
  geom_vline(xintercept = 1982,color="green", linetype="dotted") +
  geom_vline(xintercept = 2630,color="black", linetype="dotted") +
  geom_vline(xintercept = 2677,color="black", linetype="dotted")
