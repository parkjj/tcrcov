library("Rtsne")
library("rsvd")
library("data.table")
library("ggplot2")
library("Rmagic")
library("ggrepel")
library("gplots")
require("lattice")
library("dplyr")
library("Seurat")
library("patchwork")
library("pheatmap")
library("viridis")
library("scales")
library("RColorBrewer")
library("R.utils")
library("colorRamps")
library("wesanderson")
set.seed(20)


meta <- fread("../data/E-MTAB-9357.sdrf.txt")
meta <- meta[-which(meta$`Source Name` == "INCOV145-AC"),]
files <- meta[,29]
files <- files$`Derived Array Data File`
sampid <- meta$`Source Name`

files_cd4_tcr <- meta[,41]
files_cd4_tcr <- files_cd4_tcr$`Derived Array Data File`
#sampid <- sampid[1:5]

data_tcr <- NULL
SO_objects <- NULL
for(i in 1:length(sampid)){
#for(i in 265:269){
  filename_tcr <- files_cd4_tcr[i]
  filename <- files[i]
  sampname <- sampid[i]
  print(sampname)
  data_tcr_i <- fread(paste0("../data/E-MTAB-9357.processed.3/", filename_tcr))
  data_tcr_i$sample <-sampname
  data_tcr_i <- data_tcr_i[-which(data_tcr_i$TRB_1_cdr3 == "None"),]
  data_tcr <- rbind(data_tcr, data_tcr_i)
  
  data <- fread(paste0("data/pooled_gex/", filename))
  data <- as.data.frame(data)
  rownames(data) <- data$V1
  data <- data[,-1]
  data <- t(data)
  data <- exp(data)-1
  data <- data[,which(colnames(data) %in% data_tcr_i$V1)]
  
  object <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, project = sampname)
  assign(paste0(sampname), object)
  SO_objects <- c(SO_objects, paste0(sampname))
}

SO.big <- merge(get(sampid[1]), y = sapply(sampid[-1], get), add.cell.ids = SO_objects, project = "merged")
saveRDS(SO.big, file =paste0("./merged_SO_cd4.rds"))



