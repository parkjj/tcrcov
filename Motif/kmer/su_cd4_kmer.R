library("data.table")
library("ggplot2")
library("umap")
library("Rtsne")
library("locStra")
library("Matrix")
library("scales")
library("lubridate")
library("stringr")
library("RColorBrewer")
library("dplyr")
library("plyr")
library("immunarch")
library("stringi")
library("plot3D")
library("umap")
library("factoextra")
library("ggrepel")
library("pheatmap")
set.seed(20)


#################### 
# CD4 analysis
#################### 

# Note: segments of code may be rerun with different k-mer lengths

metadata <- fread("../vdj_datasets/su_cd4_reform/metadata.txt")

file.names <- dir("../vdj_datasets/su_cd4_reform", pattern =".txt")
file.names <- file.names[-which(file.names == "metadata.txt" )]

kmersize <- 7
kmermaster <- data.frame(Kmer = c(NA))
for (i in 1:length(file.names)){
  print(i)
  file <- fread(paste0("../vdj_datasets/su_cd4_reform/", file.names[i]))
  sampleid <- substr(file.names[i],1,nchar(file.names[i])-4)
  fileseqs <- NULL
  for(i in 1:dim(file)[1]){
    fileseqs <- c(fileseqs, rep(file$CDR3aa[i], file$count[i]))
  }
  seqsdf <- data.frame(CDR3.aa  = fileseqs)
  seqsdf$CDR3.aa <- as.character(seqsdf$CDR3.aa)
  seqsdf <- as_tibble(seqsdf)
  seqsdfkmers <- getKmers(seqsdf, kmersize)
  colnames(seqsdfkmers)[2] <- sampleid
  kmermaster <- merge(x = kmermaster, y = seqsdfkmers, by = "Kmer", all = TRUE)
}
kmermaster <- kmermaster[-which(is.na(kmermaster$Kmer)),]
kmermaster[is.na(kmermaster)] <- 0

fwrite(kmermaster,"su_cd4_7mers.csv", quote = F)

#### PCA ####
kmersize <- 4
kmermaster <- fread(paste0("su_cd4_",kmersize,"mers.csv"))
varexp <- apply(as.data.frame(kmermaster[,-1]), 1, var)
kmermaster <- kmermaster[order(varexp, decreasing = T),]

pcinput <-kmermaster[,-1]
pcinput <- t(t(pcinput)/colSums(pcinput))
pcinput <- as.data.frame(t(pcinput))

pcoutput <- prcomp(pcinput,center = TRUE)
summary(pcoutput)
pc <- pcoutput$x
pc <- data.frame(rownames(pc), pc[,1], pc[,2])
colnames(pc) <- c("sample", "PC1", "PC2")
pc <- cbind(pc, d_con = metadata$`Who Ordinal Scale`[match(pc$sample, metadata$Sample)])


pc$d_con[which(pc$d_con == "")] <- "HD"
pc$d_con[union(union(which(pc$d_con == "1"), which(pc$d_con == "1 or 2")), which(pc$d_con == "2"))] <- "Mild"
pc$d_con[union(which(pc$d_con == "3"), which(pc$d_con == "4"))] <- "Moderate"
pc$d_con[union(union(which(pc$d_con == "5"), which(pc$d_con == "6")), which(pc$d_con == "7"))] <- "Severe"

write.table(pc, paste0("su_cd4_pca_",kmersize,"mer.csv"), quote = F, sep = ",", row.names = F, col.names = T)

qcpc <- ggplot(data=pc, aes(x=PC1, y=PC2)) + geom_point( aes(fill=d_con), color="black", size = 5,pch=21, stroke = 1) # + scale_fill_manual(values=c("#9399ff", "#a9fffd", "#46cdcf",  "#f67280", "#c06c84", "#6c5b7b", "#355c7d"))# + #+ ggtitle("Boxplot of the sgRNA representation")
qcpc <- qcpc + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size = 0.8), axis.text=element_text(size=15), axis.title=element_text(size=20), axis.line = element_line(size = 0),  legend.title=element_text(size=0), legend.text=element_text(size=20), plot.title = element_text(size=22))
qcpc <- qcpc + labs(x = "PC1", y = "PC2") 
qcpc <- qcpc + scale_shape_manual(values=c(21, 2, 5)) + ggtitle(paste0("PCA ISB-S CD4 dataset, ",kmersize,"-mers"))
pdf(paste0("su_cd4_pca_",kmersize,"mer.pdf"), height = 5, width = 8)
qcpc
dev.off()
pdf(paste0("su_cd4_pca_",kmersize,"mer_scree.pdf"), height = 5, width = 7)
fviz_eig(pcoutput, addlabels = T, barfill = "#6E85B2",barcolor = "#261C2C") + theme(text = element_text(size = 16), axis.title = element_text(size = 16), axis.text = element_text(size = 16))+ ggtitle(paste0("Scree plot, ISB-S CD4 dataset, ",kmersize,"-mers"))
dev.off()




#### Heatmap ###########
wtsamples <- metadata$Sample[which(metadata$`Who Ordinal Scale` == "")]
expsamples1 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "1"), which(metadata$`Who Ordinal Scale` == "1 or 2")), which(metadata$`Who Ordinal Scale` == "2"))]
expsamples2 <- metadata$Sample[union(which(metadata$`Who Ordinal Scale` == "3"), which(metadata$`Who Ordinal Scale` == "4"))]
expsamples3 <- metadata$Sample[union(union(which(metadata$`Who Ordinal Scale` == "5"), which(metadata$`Who Ordinal Scale` == "6")), which(metadata$`Who Ordinal Scale` == "7"))]
expsamples4 <- union(union(expsamples1, expsamples2), expsamples3)
expid1 <- "mild"
expid2 <- "moderate"
expid3 <- "severe"
expid4 <- "pooled"
expsamples <-expsamples1
expid <- expid1

set.seed(21)
kmersize <- 3
kmermaster <- fread(paste0("su_cd4_",kmersize,"mers.csv"))
kmermaster <- as.data.frame(kmermaster)
wtsamples <- sample(wtsamples, 16)
expsamples1 <- sample(expsamples1, 16)
expsamples2 <- sample(expsamples2, 16)
expsamples3 <- sample(expsamples3, 16)
expsamples4 <- union(union(expsamples1, expsamples2), expsamples3)
kmermaster <- cbind(kmermaster[,1], 
                    kmermaster[,which(colnames(kmermaster) %in% wtsamples)],
                    kmermaster[,which(colnames(kmermaster) %in% expsamples1)],
                    kmermaster[,which(colnames(kmermaster) %in% expsamples2)],
                    kmermaster[,which(colnames(kmermaster) %in% expsamples3)])

varexp <- apply(as.data.frame(kmermaster[,-1]), 1, var)
kmermaster <- kmermaster[order(varexp, decreasing = T),]
kmermaster <- kmermaster[c(1:40),]
kmernorm <-kmermaster[,-1]
kmernorm <- t(t(kmernorm)/colSums(kmernorm))
kmernorm <- as.data.frame(kmernorm)

zscore <- as.data.frame(kmernorm)
zscore <- log2(zscore+1)
zscore <- t(scale(t(zscore)))

zscore.plt <- zscore#[match(goi, rownames(zscore)),]
zscore.plt <- data.frame(kmer = kmermaster[,1],zscore.plt)
rownames(zscore.plt) <- NULL
roworderindex1 <- rowMeans(zscore.plt[,which(colnames(zscore.plt) %in% wtsamples)])
roworderindex2 <- rowMeans(zscore.plt[,which(colnames(zscore.plt) %in% expsamples4)])
roworderindex3 <- roworderindex2-roworderindex1
zscore.plt <- zscore.plt[order(roworderindex3, decreasing = F),]
zscore.plt$kmer <- factor(zscore.plt$kmer, levels = zscore.plt$kmer)
phplt <- as.matrix(zscore.plt[,-1])
rownames(phplt) <- zscore.plt$kmer 
pdf(paste0("dex/hm_ph_su_cd4_",kmersize,"mer.pdf"), width = 16, height = 12)
pheatmap(phplt, cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(100), border_color = "grey90")
dev.off()
zscore.plt <- melt(zscore.plt)
zscore.plt$variable <- factor(as.character(zscore.plt$variable), levels = c(wtsamples, expsamples1, expsamples2, expsamples3))
plt <- ggplot(zscore.plt, aes(variable, kmer)) + geom_tile(aes(fill = value)) +scale_fill_gradientn(colours=rev(brewer.pal(n = 10, name = "RdBu")))#+scale_fill_viridis(option="inferno")#+scale_fill_gradientn(colours=rev(brewer.pal(n = 10, name = "RdBu")))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=12), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=10)) 
plt <- plt+ labs(x = "", y = "")
plt <- plt + scale_y_discrete(position = "right")
pdf(paste0("dex/hm_su_cd4_",kmersize,"mer.pdf"), width = 16, height = 12)
plt
dev.off()
