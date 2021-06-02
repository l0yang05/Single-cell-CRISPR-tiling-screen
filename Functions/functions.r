## helper functions
mk.Seurat <- function(dataDir, projectID,genome="human",min.cells=5,min.genes = 200,runVinPlot= T,mksub = F, nFeature.low =NULL,nFeature.high=NULL, percentMT.cutoff=NULL,lib.type="Expression only"){
  # create seurat object and make QC plots
  dat <- Read10X(data.dir=dataDir)
  print("read in 10x data.......")
  if(lib.type=="Expression only"){
    print("Provided 10x dataset is based on EXpression only library")
    dat.seurat <- CreateSeuratObject(counts = dat, min.cells = min.cells, min.features = min.genes, project = projectID)

  }
  else if(lib.type=="Expression plus CRISPR"){
    print("Provided 10x dataset is based on feature bacrode with EXpression and CRISPR both!")
    dat.exp <- dat$`Gene Expression`
    dat.seurat <- CreateSeuratObject(counts = dat.exp, min.cells = min.cells, min.features = min.genes, project = projectID)

  }
  #dat.seurat <- CreateSeuratObject(counts = dat, min.cells = min.cells, min.features = min.genes, project = projectID)
  print("Seurat object has been created.......")
  ## add mt% to meta.data column
  if(genome=="human"){
    dat.seurat[["percent.mt"]] <- PercentageFeatureSet(dat.seurat, pattern = "^MT-")
  }
  else if(genome=="mouse"){dat.seurat[["percent.mt"]] <- PercentageFeatureSet(dat.seurat, pattern = "^mt-")}

  p1 <- VlnPlot(dat.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(runVinPlot)
  plotName <- paste("VinPlot_",projectID,".pdf",sep="")
  print(plotName)
  if(runVinPlot){
    pdf(plotName, width = 12, height = 8)
    print(p1)
    dev.off()
  }

  if (mksub){
    print("start subsetting dataset")
    nFeature.low <<- nFeature.low  ## reassign global variable
    nFeature.high <<- nFeature.high
    percentMT.cutoff <<- percentMT.cutoff
    dat.sub <- subset(dat.seurat, subset = nFeature_RNA > nFeature.low & nFeature_RNA < nFeature.high & percent.mt < percentMT.cutoff)
    nCells <- nrow(dat.sub@meta.data)
    fname <- paste(projectID,"_",nCells,".rd",sep="")
    print(fname)
    saveRDS(dat.sub,file=fname)


  }

}

add_category <-function(feature_call, anno.tbl){
  ## if feature_call does not exisit in the anno,tbl, add feature_call as a new category
  idx = match(feature_call, anno.tbl$id)
  if(is.na(idx)){
    cat.info<-feature_call
  }
  else{cat.info <- anno.tbl$category[idx]}
  cat.info
}

##help functions of cor-mat analysis
## input count matrix(data.frame): 1st column is "pheno.feature_call" or other ids to be grouped by,rownames are cell id
## output ， matrix, columns are guides ordered by position, rows are selected genes
ctrlIDs <- c("LUC_474", "Ren_55", "Rosa26", "tRFP657_294")
cormat_analysis <- function(mat, value_to_cor = "median", ctrlIDs){
  ##1 generate matrix to be used for plotting heatmap, either median or mean
  group_by_id <- colnames(mat)[1] # eg. pheno.feature_call
  mat.v1 <- mat %>% group_by_(.dots = group_by_id) %>%
      summarise_all(value_to_cor)
  ##2 add pos id based on group_by_id: eg. pheno.feature_call
  mutName = "pos"   ## add position info
  mutf=paste0('as.numeric(word(',group_by_id,',-2,sep=fixed(\'_\')))')
  mat.v1.pos <- mat.v1 %>% mutate_(.dots = setNames(mutf, mutName)) %>%
    arrange(pos)
  ##4 remove ctrl IDs, mat.v1.pos is tibble
  mat.v2.pos <- mat.v1.pos %>%
    filter(!mat.v1.pos[[1]]%in%ctrlIDs)
  rnames <- as.character(mat.v2.pos[[1]])
  ## remove first column and last position column
  mat.v3.pos <- mat.v2.pos[,-c(1,ncol(mat.v2.pos))]
  row.names(mat.v3.pos) <- rnames
  mat.v3.pos.tr <- t(mat.v3.pos)


  # Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

## smoothen corr ratio over position

smoothen.cor <- function(rowdat,pos,bw=NULL){
  ## input is two-column tbl from one row of the cor-matrix, pos and cor
  # length is how many data points to have
  dat <- data.frame(pos = pos, cor = rowdat)
  if(is.null(bw)){
    bw=round(max(diff(pos)))
  }
  n = dat$pos[length(pos)]-dat$pos[1] + 1
  print(n)
  dat.sm <- ksmooth(dat$pos, dat$cor, kernel="normal",bandwidth=0.5*bw,n.points = n)
  dat.sm.all <- data.frame(pos=round(dat.sm$x), cor=dat.sm$y)
  temp <-dat.sm.all$cor
  names(temp) <- dat.sm.all$pos
  temp

}
