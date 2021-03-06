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
library("VennDiagram")
library("ggsignif")
library("ggrepel")

set.seed(20)


#################### 
# Adaptive vgene
#################### 

metaHD <- fread("../vdj_datasets/adaptiveHD/metadata.txt")
vgenepropHD <- fread("vgeneprop_vdj_adaptiveHD.csv")
vgenepropHD <- as.data.frame(vgenepropHD)

meta <- fread("../vdj_datasets/adaptive/metadata.txt")
vgeneprop <- fread("vgeneprop_vdj_adaptive.csv")
vgeneprop <- as.data.frame(vgeneprop)


# mean prop
removegenes <- c("TRBV12-3/12-4", "TRBV20-or09_02", "TRBV23-or09_02", "TRBV25-or09_02", "TRBV26-or09_02", "TRBV3-1/03-2", "TRBV6-2/06-3", "TRBVA", "TRBVA-1",
                 "TRBV21-or09_02", "TRBV22-or09_02", "TRBV24-or09_02", "TRBV29-or09_02", "TRBVA-or09_02", "unknown")

vgenepropHD <- vgenepropHD[-which(vgenepropHD$vgene %in% removegenes),]
vgeneprop <- vgeneprop[-which(vgeneprop$vgene %in% removegenes),]

reindexHD <- substr(vgenepropHD$vgene, 5, nchar(vgenepropHD$vgene))
reindexHD <- as.numeric(gsub("-",".",reindexHD))
reindexHD <- order(reindexHD)

reindex <- substr(vgeneprop$vgene, 5, nchar(vgeneprop$vgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)

#vgenepropHD_vgene <- vgenepropHD$vgene[reindexHD]
#vgeneprop_vgene <- vgeneprop$vgene[reindex]

vgeneprop_HD <- vgenepropHD[reindexHD,]
vgeneprop_COV <- vgeneprop[reindex,]

vgeneprop_HD[,-1] <- t(t(vgeneprop_HD[,-1]) / colSums(vgeneprop_HD[,-1]))
vgeneprop_COV[,-1] <- t(t(vgeneprop_COV[,-1]) / colSums(vgeneprop_COV[,-1]))
write.csv(vgeneprop_HD,"vgeneprop_vdj_adaptiveHD_nl.csv" ,quote = F,  row.names = T)
write.csv(vgeneprop_COV,"vgeneprop_vdj_adaptive_nl.csv" ,quote = F,  row.names = T)



pltd_HD_melt <- melt(cbind(rep("HD",dim(vgeneprop_HD)[1]), vgeneprop_HD))
pltd_COV_melt <- melt(cbind(rep("COVID",dim(vgeneprop_COV)[1]), vgeneprop_COV))

colnames(pltd_HD_melt) <- c("condition", "vgene", "sample", "prop")
colnames(pltd_COV_melt) <- c( "condition", "vgene","sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_COV_melt)
pltd$vgene <- factor(pltd$vgene, levels = unique(pltd$vgene))
pltd$condition <- factor(pltd$condition, levels = unique(pltd$condition))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=0.2, alpha=0.5)+ scale_color_manual(values = c("#787A91", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white",lwd=0.3, outlier.colour=NA, position=position_dodge(width=0.9))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits=c(0,60))
plt <- plt + labs(x= "", y = "Percentage")
pdf("adaptive_vgene.pdf", height = 6.5, width = 15)
plt
dev.off()
png("adaptive_vgene.png", height = 6.5, width = 15,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=vgene, y=prop, color=condition)) + geom_point(pch=26, position=position_jitterdodge(dodge.width=0.9), size=0.2, alpha=0.5)+ scale_color_manual(values = c("#787A91", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits=c(0,60))
plt <- plt + labs(x= "", y = "Percentage")
pdf("adaptive_vgene.fo.pdf", height = 6.5, width = 15)
plt
dev.off()


# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$vgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,61)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 61)[1]]
  part_covid <- parts[seq(i, length(parts), 61)[2]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_covid)])
  
  pvalouts <- c(wout1$p.value)
  print(pvalouts)
  pvalues <- c(pvalues, pvalouts)
}

pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))
statvals <- as.data.frame(pvalues)
rownames(statvals) <- as.character(unique(pltd$vgene))
write.csv(statvals,"adaptive_vgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"adaptive_vgene_wilcoxon_pvalues_sig_5e-2.csv" ,quote = F,  row.names = T)
write.csv((statvals < 1e-4)*1,"adaptive_vgene_wilcoxon_pvalues_sig_1e-4.csv" ,quote = F,  row.names = T)
write.csv((statvals < 1e-6)*1,"adaptive_vgene_wilcoxon_pvalues_sig_1e-6.csv" ,quote = F,  row.names = T)



vgeneset <- vgeneprop_HD$vgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(vgeneset)){
  print(i)
  set1 <- as.numeric(vgeneprop_COV[i,-1])
  set2 <- as.numeric(vgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_vgene <- data.frame(vgene = vgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_vgene <- allprop_stats_vgene[order(allprop_stats_vgene$sig, decreasing = T),]
write.csv(allprop_stats_vgene,"allprop_stats_vgene.csv" ,quote = F,  row.names = F)



## Adaptive jgene
metaHD <- fread("../vdj_datasets/adaptiveHD/metadata.txt")
jgenepropHD <- fread("jgeneprop_vdj_adaptiveHD.csv")
jgenepropHD <- as.data.frame(jgenepropHD)

meta <- fread("../vdj_datasets/adaptive/metadata.txt")
jgeneprop <- fread("jgeneprop_vdj_adaptive.csv")
jgeneprop <- as.data.frame(jgeneprop)


removegenes <- c("TRBJ2-2P")

jgenepropHD <- jgenepropHD[-which(jgenepropHD$jgene %in% removegenes),]
jgeneprop <- jgeneprop[-which(jgeneprop$jgene %in% removegenes),]

reindexHD <- substr(jgenepropHD$jgene, 5, nchar(jgenepropHD$jgene))
reindexHD <- as.numeric(gsub("-",".",reindexHD))
reindexHD <- order(reindexHD)

reindex <- substr(jgeneprop$jgene, 5, nchar(jgeneprop$jgene))
reindex <- as.numeric(gsub("-",".",reindex))
reindex <- order(reindex)


jgeneprop_HD <- jgenepropHD[reindexHD,]
jgeneprop_COV <- jgeneprop[reindex,]

jgeneprop_HD[,-1] <- t(t(jgeneprop_HD[,-1]) / colSums(jgeneprop_HD[,-1]))
jgeneprop_COV[,-1] <- t(t(jgeneprop_COV[,-1]) / colSums(jgeneprop_COV[,-1]))
write.csv(jgeneprop_HD,"jgeneprop_vdj_adaptiveHD_nl.csv" ,quote = F,  row.names = T)
write.csv(jgeneprop_COV,"jgeneprop_vdj_adaptive_nl.csv" ,quote = F,  row.names = T)




pltd_HD_melt <- melt(cbind(rep("HD",dim(jgeneprop_HD)[1]), jgeneprop_HD))
pltd_COV_melt <- melt(cbind(rep("COVID",dim(jgeneprop_COV)[1]), jgeneprop_COV))

colnames(pltd_HD_melt) <- c("condition", "jgene", "sample", "prop")
colnames(pltd_COV_melt) <- c( "condition", "jgene","sample", "prop")
pltd <- rbind(pltd_HD_melt, pltd_COV_melt)
pltd$jgene <- factor(pltd$jgene, levels = unique(pltd$jgene))
pltd$condition <- factor(pltd$condition, levels = unique(pltd$condition))
pltd$prop <- pltd$prop*100

plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(position=position_jitterdodge(dodge.width=0.9), size=0.2, alpha=0.5)+ scale_color_manual(values = c("#787A91", "#CD113B"))
plt <- plt + geom_boxplot(fill = "white", outlier.colour=NA, position=position_dodge(width=0.9))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits=c(0,70))
plt <- plt + labs(x= "", y = "Percentage")
pdf("adaptive_jgene.pdf", height = 6, width = 9)
plt
dev.off()
png("adaptive_jgene.png", height = 6, width = 9,units = "in", res = 300)
plt
dev.off()

plt <- ggplot(pltd, aes(x=jgene, y=prop, color=condition)) + geom_point(pch=26,position=position_jitterdodge(dodge.width=0.9), size=0.2, alpha=0.5)+ scale_color_manual(values = c("#787A91", "#CD113B"))
plt <- plt + theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=2),axis.line = element_blank(), axis.text.x = element_text(size=20,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=20), axis.title=element_text(size=25), plot.title=element_text(size=30),legend.position = "none", legend.title=element_text(size=0), legend.text=element_text(size=8)) 
plt <- plt + scale_y_continuous(breaks=pretty_breaks(6), limits=c(0,70))
plt <- plt + labs(x= "", y = "Percentage")
pdf("adaptive_jgene.fo.pdf", height = 6, width = 9)
plt
dev.off()



# stats
pvalues <- NULL
pltd$partition <- paste0(pltd$jgene, "_", pltd$condition)
parts <- as.character(unique(pltd$partition))

for(i in seq(1,13)){
  print(i)
  part_hd <- parts[seq(i, length(parts), 13)[1]]
  part_covid <- parts[seq(i, length(parts), 13)[2]]
  
  wout1 <- wilcox.test(pltd$prop[which(pltd$partition == part_hd)],
                       pltd$prop[which(pltd$partition == part_covid)])
  
  pvalouts <- c(wout1$p.value)
  print(pvalouts)
  pvalues <- c(pvalues, pvalouts)
}

pvalues <- p.adjust(pvalues, method = "BH", n = length(pvalues))

statvals <- as.data.frame(pvalues)
rownames(statvals) <- as.character(unique(pltd$jgene))
write.csv(statvals,"adaptive_jgene_wilcoxon_pvalues.csv" ,quote = F,  row.names = T)
write.csv((statvals < 0.05)*1,"adaptive_jgene_wilcoxon_pvalues_sig_5e-2.csv" ,quote = F,  row.names = T)
write.csv((statvals < 1e-4)*1,"adaptive_jgene_wilcoxon_pvalues_sig_1e-4.csv" ,quote = F,  row.names = T)
write.csv((statvals < 1e-6)*1,"adaptive_jgene_wilcoxon_pvalues_sig_1e-6.csv" ,quote = F,  row.names = T)



jgeneset <- jgeneprop_HD$jgene
w <- NULL
m1 <- NULL
m2 <- NULL
for(i in 1:length(jgeneset)){
  print(i)
  set1 <- as.numeric(jgeneprop_COV[i,-1])
  set2 <- as.numeric(jgeneprop_HD[i,-1])
  wout <- wilcox.test(set1, set2)
  w <- c(w, wout$p.value)
  m1 <- c(m1, mean(set1))
  m2 <- c(m2, mean(set2))
}
w_adj <- p.adjust(w, "BH")
w_adj_sig <- -log10(w_adj)
md <- m1-m2

allprop_stats_jgene <- data.frame(jgene = jgeneset,
                                  mean_covid = m1,
                                  mean_hd = m2,
                                  mean_diff = md,
                                  logfc = log(m1/m2),
                                  pval = w,
                                  padj = w_adj,
                                  sig = w_adj_sig)
allprop_stats_jgene <- allprop_stats_jgene[order(allprop_stats_jgene$sig, decreasing = T),]
write.csv(allprop_stats_jgene,"allprop_stats_jgene.csv" ,quote = F,  row.names = F)

