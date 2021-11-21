library("Rtsne")
library("rsvd")
library("data.table")
library("ggplot2")
#library("Rmagic")
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
set.seed(20)


i <-"merged"
load(paste0(i, "_workspace.RData"))


object <- SO.big
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 1000000)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(object), 10)
plot1 <- VariableFeaturePlot(object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))

object <- FindNeighbors(object, dims = 1:10)
object <- FindClusters(object, resolution = 0.6)
object <- RunUMAP(object, dims = 1:10)

object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(object.markers, paste0("dex_allmarkers.csv"), quote = F)

top10 <- object.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

out <- Embeddings(object = object, reduction = "pca")
fwrite(out, "pca_embeddings_corrected.csv")
out <- Embeddings(object = object, reduction = "umap")
fwrite(out, "umap_embeddings_corrected.csv")
out <- object@meta.data
fwrite(out, "meta_corrected.csv", row.names = T)
#out <- GetAssayData(object = object, slot = "data")
#out <- as.data.frame(as.matrix(out))
#fwrite(out, "logcounts.csv", row.names = T)


####################################################################
############# Save and load workspace ############# 
####################################################################
quickstats <- dim(object)
quickstats <- c(quickstats, length(levels(object)))
names(quickstats) <- c("features", "samples", "clusters")
write.table(quickstats, paste0("quickstats.csv"), quote = F, col.names = F, sep = ",")
write.table(levels(object), paste0("label_template.csv"), quote = F, col.names = F,  row.names = F, sep = ",")
write.table(all.genes, paste0("allgenes.csv"), quote = F, col.names = F,  row.names = F, sep = ",")

save.image(file =paste0(i, "_workspace_reprocessed.RData"))