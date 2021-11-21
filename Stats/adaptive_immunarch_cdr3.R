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
library("ggrepel")
library("car")
library("lmPerm")
library("VennDiagram")

set.seed(20)


#################### 
# CD4 analysis total clonotypes
#################### 


meta <- fread("../vdj_datasets/adaptive/metadata.txt")
metaHD <- fread("../vdj_datasets/adaptiveHD/metadata.txt")

allpropCOV <- fread("allprop_cdr3_adaptive.csv")
allpropCOV <- as.data.frame(allpropCOV)
allpropCOV_hold <- allpropCOV
allpropHD <- fread("allprop_cdr3_adaptiveHD.csv")
allpropHD <- as.data.frame(allpropHD)
allpropHD_hold <- allpropHD

#allprop_meanprop <- rowMeans(allprop[,-1])
#allpropHD_meanprop <- rowMeans(allpropHD[,-1])

allpropCOV <- allpropCOV_hold
allpropHD <-allpropHD_hold
allpropCOV <-  allpropCOV[which(rowSums(allpropCOV[,-1]) > 0.005),]
allpropHD <-  allpropHD[which(rowSums(allpropHD[,-1]) > 0.005),]

cdr3set <- union(allpropHD$cdr3, allpropCOV$cdr3)
allprop <- merge(x = allpropHD, y = allpropCOV, by = "cdr3", all = TRUE)
allprop[is.na(allprop)] <- 0

allpropCOV <-allprop[,c(1,90:1564)]
allpropHD <- allprop[,1:89]
cdr3set <- allpropHD$cdr3


q <- data.frame(cdr3 = allpropHD$cdr3 ,
                diff = (rowMeans(allpropCOV[,-1]) - rowMeans(allpropHD[,-1])))
q <- q[order(q$diff),]
q$cdr3 <- factor(as.character(q$cdr3), levels = as.character(q$cdr3))
q$hit <- rep("CDR3", dim(q)[1])
q$hit[1:7] <- "HD"
q$hit[(length(q$hit)-6):length(q$hit)] <- "COVID-19"

plt <- ggplot(q, aes(x=cdr3, y=diff, color=hit)) + geom_point( size = 2) + scale_colour_manual(values = c("gray","#A239EA","#11cbd7"))# + scale_colour_manual(values = c("gray","#11cbd7","#7e6bc4",  "#fa4659"))
plt <- plt + theme_classic() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_text(size=25),axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=30), axis.line = element_line(colour = "black", size = 0.8), plot.title=element_text(size=30), legend.title=element_text(size=0), legend.text=element_text(size=15))
plt <- plt + labs(x = "CDR3", y = expression(Delta*" Mean Proportion"))
plt <- plt + geom_point(data=subset(q, hit == "HD"), size = 5, show.legend = F) 
plt <- plt + geom_point(data=subset(q, hit == "COVID-19"), size = 5, show.legend = F) 
plt <- plt + geom_text_repel(data=subset(q, hit == "HD"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3), size = 7,show.legend = FALSE, nudge_x = 1000, hjust = 0, direction = "y")
plt <- plt + geom_text_repel(data=subset(q, hit == "COVID-19"), mapping=aes(x=cdr3, y=diff, color=hit, label=cdr3),nudge_y= 5e-4, size = 7,show.legend = FALSE, nudge_x = -1000, hjust = 1, direction = "y")
plt <- plt + scale_y_continuous(breaks = pretty_breaks(5), limits = c(-2.6E-3, 2e-3))
pdf("waterfall_covidvshd.pdf", height = 6, width = 10)
plt
dev.off()





allpropCOV <- allpropCOV_hold
allpropHD <-allpropHD_hold

cdr3set <- union(allpropHD$cdr3, allpropCOV$cdr3)
allprop <- merge(x = allpropHD, y = allpropCOV, by = "cdr3", all = TRUE)
allprop[is.na(allprop)] <- 0

allpropCOV <-allprop[,c(1,90:1564)]
allpropHD <- allprop[,1:89]
cdr3set <- allpropHD$cdr3



ns <- 15
pltd_HD <- data.frame(cdr3 = allpropHD$cdr3, meanprop = rowMeans(allpropHD[, -1]))
pltd_HD$meanprop <- pltd_HD$meanprop*100
pltd_HD <- pltd_HD[order(pltd_HD$meanprop, decreasing = T),]
pltd_HD <- pltd_HD[1:ns,]
pltd_HD$group <- rep("HD", ns)
pltd_HD$cdr3 <- factor(pltd_HD$cdr3, levels = as.character(pltd_HD$cdr3))
plt <- ggplot(pltd_HD, aes(x=cdr3, y=meanprop)) + geom_bar(stat = "identity", color="black", fill="#B2B1B9") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("Healthy Donor")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_mean_hd.pdf", height = 8, width = 6)
plt
dev.off()


pltd_COV <- data.frame(cdr3 = allpropCOV$cdr3, meanprop = rowMeans(allpropCOV[, -1]))
pltd_COV$meanprop <- pltd_COV$meanprop*100
pltd_COV <- pltd_COV[order(pltd_COV$meanprop, decreasing = T),]
pltd_COV <- pltd_COV[1:ns,]
pltd_COV$group <- rep("COVID-19", ns)
pltd_COV$cdr3 <- factor(pltd_COV$cdr3, levels = as.character(pltd_COV$cdr3))
plt <- ggplot(pltd_COV, aes(x=cdr3, y=meanprop)) + geom_bar(stat = "identity", color="black", fill="#A239EA") 
plt <- plt + theme_bw() + theme(panel.background = element_rect(colour = "black", size=2, fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=20, colour = "black"), axis.text.x=element_text(size=20, colour = "black", angle = 90, vjust = 0.5, hjust=1), axis.title=element_text(size=20), axis.line = element_line(colour = "black", size = 0), plot.title=element_text(size=20), legend.title=element_text(size=0), legend.text=element_text(size=15)) 
plt <- plt + labs(x = "", y = "Mean CDR3 usage (%)")
plt <- plt + ggtitle("COVD-19")# + guides(colour = guide_legend(override.aes = list(size=5)))
plt <- plt + ylim(0,0.5) + coord_fixed(ratio = 20)
pdf("cdr3_mean_covid.pdf", height = 8, width = 6)
plt
dev.off()


allprop_HD <- fread("allprop_cdr3_su_cd4_processed_HD.csv")
allprop_mild <- fread("allprop_cdr3_su_cd4_processed_mild.csv")
allprop_mod <- fread("allprop_cdr3_su_cd4_processed_mod.csv")
allprop_severe <- fread("allprop_cdr3_su_cd4_processed_severe.csv")
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.0001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.0001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.0001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.0001)]))
pattern_pool_cd4 <- union(union(pattern_mild, pattern_mod), pattern_sev)


allprop_HD <- fread("allprop_cdr3_su_cd8_processed_HD.csv")
allprop_mild <- fread("allprop_cdr3_su_cd8_processed_mild.csv")
allprop_mod <- fread("allprop_cdr3_su_cd8_processed_mod.csv")
allprop_severe <- fread("allprop_cdr3_su_cd8_processed_severe.csv")
pattern_hd <- unique(as.character(allprop_HD$cdr3[which(allprop_HD$meanprop > 0.0001)]))
pattern_mild <- unique(as.character(allprop_mild$cdr3[which(allprop_mild$meanprop > 0.0001)]))
pattern_mod <- unique(as.character(allprop_mod$cdr3[which(allprop_mod$meanprop > 0.0001)]))
pattern_sev <- unique(as.character(allprop_severe$cdr3[which(allprop_severe$meanprop > 0.0001)]))
pattern_pool_cd8 <- union(union(pattern_mild, pattern_mod), pattern_sev)


#allpropCOV_means <- rowMeans(allpropCOV[, -1]) 
pattern_covid <- unique(as.character(allpropCOV$cdr3[which(allpropCOV_means > 0.00001)]))


pdf("intersect_analysis.pdf", height = 6, width = 6)
draw.triple.venn(area1 = length(pattern_pool_cd4), area2 = length(pattern_pool_cd8), area3 = length(pattern_covid), 
                 n12 = length(intersect(pattern_pool_cd4, pattern_pool_cd8)), 
                 n23 = length(intersect(pattern_pool_cd8, pattern_covid)), 
                 n13 = length(intersect(pattern_pool_cd4, pattern_covid)), 
                 n123 = length(intersect(intersect(pattern_pool_cd4, pattern_pool_cd8), pattern_covid)), 
                 cex = 2,
                 #cat.cex = 2,
                 category = c("ISB-S CD4", "ISB-S CD8", "AB COVID-19"),lwd = 3, col = "black",
                 fill = c("#FFE3FE", "#B4AEE8", "#93329E"));
dev.off()




