#!/usr/bin/Rscript --vanilla

####import libraries and data####
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(scDblFinder)

#####import data (source code deleted/ edited for privacy from the original script)#####

# Set path to data directory

# Create a seurat object for each regions

#####quality control#####

# Merge the data  
merged_data <- merge(Amygdala, y = c(Cerebellum, Hippocampus, Occipital, Parietal),
                     add.cell.ids = c("Amygdala", "Cerebellum", "Hippocampus", "Occipital", "Parietal"),
                     project = "")

merged_data$sample　<- rownames(merged_data@meta.data)
merged_data@meta.data <- separate(merged_data@meta.data, col = "sample", 
                                  into = c("region", "Barcode"), sep = "_")

#check to make sure if it has merged correctly
#View(merged_data@meta.data)
#unique(merged_data@meta.data$region)


#create a column for mtRNA
merged_data[["percent.mt"]] <- PercentageFeatureSet(merged_data, pattern = "^MT-")


vplotQC <- VlnPlot(object = merged_data,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  group.by = "region",
                  pt.size = 0)


# Apply QC threshold
merged_data_filtered <- subset(merged_data, subset = nFeature_RNA > 200 & percent.mt < 5)


#Check to see if there is any batch effects
# merged_data_filtered <- NormalizeData(merged_data_filtered)
# merged_data_filtered <- FindVariableFeatures(merged_data_filtered)
# merged_data_filtered <- ScaleData(merged_data_filtered)
# merged_data_filtered <- RunPCA(merged_data_filtered)
# 
# merged_data_filtered <- FindNeighbors(merged_data_filtered, dims = 1:20)
# merged_data_filtered <- FindClusters(merged_data_filtered, resolution = 0.5)
# merged_data_filtered <- RunUMAP(merged_data_filtered, dims = 1:20)
# 
# dimplot_pre-int <- DimPlot(merged_data_filtered, reduction = "umap", group.by = "region")



#####doublet removal#####

# separate the merged data (prepare for integration)
obj.list <- SplitObject(merged_data_filtered, split.by = "region")
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
  obj.list[[i]] <- ScaleData(obj.list[[i]])
  obj.list[[i]] <- RunPCA(obj.list[[i]])
  sce <- scDblFinder(GetAssayData(obj.list[[i]], slot="counts"))
  obj.list[[i]]$scDblFinder.score <- sce$scDblFinder.score
  obj.list[[i]]$scDblFinder.class <- sce$scDblFinder.class
  obj.list[[i]] <- subset(obj.list[[i]], subset = scDblFinder.class=="singlet")
}


# violin plot to compare pre and post doublet finder
vplotC_predbl <- VlnPlot(Cerebellum, features = "nFeature_RNA", pt.size = 0, y.max = 12000) + theme(legend.position = 'none')
vplotC_postdbl <- VlnPlot(obj.list$Cerebellum, features = "nFeature_RNA", group.by = "region", pt.size = 0, y.max = 12000) + theme(legend.position = 'none')
vplotO_predbl <- VlnPlot(Occipital, features = "nFeature_RNA", pt.size = 0, y.max = 12000) + theme(legend.position = 'none')
vplotO_postdbl <- VlnPlot(obj.list$Occipital, features = "nFeature_RNA", group.by = "region", pt.size = 0, y.max = 12000) + theme(legend.position = 'none')

grid.arrange(vplotC_predbl, vplotC_postdbl, vplotO_predbl, vplotO_postdbl, ncol=4)

#vplotA1 <- VlnPlot(Amygdala, features = "nFeature_RNA")
#vplotA2 <- VlnPlot(obj.list$Amygdala, features = "nFeature_RNA", group.by = "region")
# grid.arrange(vplotA1, vplotA2, ncol=2)
#
#vplotH1 <- VlnPlot(Hippocampus, features = "nFeature_RNA")
#vplotH2 <- VlnPlot(obj.list$Hippocampus, features = "nFeature_RNA", group.by = "region")
#grid.arrange(vplotH1, vplotH2, ncol=2)
#
# vplotP1 <- VlnPlot(Parietal, features = "nFeature_RNA")
# vplotP2 <- VlnPlot(obj.list$Parietal, features = "nFeature_RNA", group.by = "region")
# grid.arrange(vplotP1, vplotP2, ncol=2)




#####Perform integration to correct for batch effects#####
# Select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

#Find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = "rpca")
to_integrate <- Reduce(intersect, lapply(anchors@object.list, rownames))

# Integrate data 
data.integrated <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data.
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.integrated) <- "integrated"



###### Run the standard workflow for visualization and clustering#####
data.integrated <- ScaleData(object = data.integrated)
data.integrated <- RunPCA(object = data.integrated)

#check dimentionality of data
ElbowPlot(data.integrated, ndims = 30)
# data.integrated <- JackStraw(data.integrated, num.replicate = 100)
# data.integrated <- ScoreJackStraw(data.integrated, dims = 1:20)
# JackStrawPlot(data.integrated, dims = 1:20)


data.integrated <- FindNeighbors(data.integrated, dims = 1:6)


######cluster resolution#####
# investigating resolution parameters
# for(i in 1:40){
#   data.integrated <- FindClusters(data.integrated, resolution = i/40)
# }

# number_of_cluster <- list()
# resolution <- list()

# for(j in 1:40){
#   if (j == 1) {
#     number_of_cluster <- append(number_of_cluster, max(as.numeric(as.character(data.integrated@meta.data[,8+j])), na.rm = TRUE)+1)
#   }
#   else{
#     number_of_cluster <- append(number_of_cluster, max(as.numeric(as.character(data.integrated@meta.data[,9+j])), na.rm = TRUE)+1) 
#   }
#    resolution <- append(resolution, j/40)
# }
# 
# plot(resolution, number_of_cluster) 


#obtains 8 clusters
#data.integrated <- FindClusters(data.integrated, resolution = 0.1)

#obtains 6 clusters
data.integrated <- FindClusters(data.integrated, resolution = 0.025)



#####UMAP visualisation#####
data.integrated <- RunUMAP(object = data.integrated, dims = 1:6)

dimplot_unlabeled <- DimPlot(data.integrated, reduction = "umap")

#QC
dimplot_region <- DimPlot(data.integrated, reduction = "umap", group.by = "region")
vplot_mt <- VlnPlot(data.integrated, features = "percent.mt", pt.size = 0)
fplot_mt <- FeaturePlot(data.integrated, features = "percent.mt")

#top 5 marker genes
data.integrated.markers <- FindAllMarkers(data.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.integrated.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  print(n=50) %>%
  capture.output(file = "FindAllMarkers.txt")


#A list of marker genes for major cell type based on current literature
exc_gene <- c("SLC17A7",	 "SATB2",	 "SNAP25",	 "FEZF2",	 "LAMP5",	 "RORB",	 "SYT1",	 "VIP")
inh_gene <- c("GAD1",	 "GAD2",	 "ADARB2",	 "LHX6",	 "PVALB",	 "SNAP25",	 "SOX6")
ast_gene <- c("AQP4",	 "GFAP",	 "FGFR3",	 "GJA1",	 "GJB6",	 "SLC1A2",	 "SLC14A1")
mic_gene <- c( "CSF1R",	 "P2RY12",	 "CD74",	"C1QA",	"C1QB",	"CTSS",	 "IRF8",	 "SPP1")
opc_gene <- c( "PDGFRA",	 "BCAN",	"CSPG4",	 "PTPRZ1",	 "VCAN")
oli_gene <- c( "MOG",	"MAG",	"PLP1",	 "CNP",	 "MBP",	 "MOBP",	 "OPALIN")
per_gene <- c( "ATP1A2",	 "CLDN5",	 "DCN",	 "EPS8",	 "FLT1",	 "ITIH5",	 "LAMA2",	 "NOTCH3")
end_gene <- c("CLDN5",	 "ACTA2",	 "EPAS1",	 "FLT1",	 "KCNJ78",	 "NOTCH3",	 "ZEB1")


#violin plots
vplot_exc <- VlnPlot(data.integrated, assay = "RNA", features = exc_gene, pt.size = 0, ncol=4)  #Excitatory Neurons
vplot_inh <- VlnPlot(data.integrated, features = inh_gene, pt.size = 0, assay = "RNA", ncol=4)  #Inhibitory Neurons
vplot_ast <- VlnPlot(data.integrated, features = ast_gene, pt.size = 0, assay = "RNA", ncol=4)  #Astrocytes
vplot_mic <- VlnPlot(data.integrated, features = mic_gene, pt.size = 0, assay = "RNA", ncol=4)  #Microglia
vplot_opc <- VlnPlot(data.integrated, features = opc_gene, pt.size = 0, assay = "RNA", ncol=4)  #Oligodendrocyte Precursor Cells (OPC)
vplot_oli <- VlnPlot(data.integrated, features = oli_gene, pt.size = 0, assay = "RNA", ncol=4)  #Oligodendrocytes
vplot_per <- VlnPlot(data.integrated, features = per_gene, pt.size = 0, assay = "RNA", ncol=4)  #Pericytes
vplot_end <- VlnPlot(data.integrated, features = end_gene, pt.size = 0, assay = "RNA", ncol=4)  #Endothelial Cells


#feature plots
fplot_exc <- FeaturePlot(data.integrated, features = exc_gene, ncol=4)  #Excitatory Neurons
fplot_inh <- FeaturePlot(data.integrated, features = inh_gene, ncol=4)  #Inhibitory Neurons
fplot_ast <- FeaturePlot(data.integrated, features = ast_gene, ncol=4)  #Astrocytes
fplot_mic <- FeaturePlot(data.integrated, features = mic_gene, ncol=4)  #Microglia
fplot_opc <- FeaturePlot(data.integrated, features = opc_gene, ncol=4)  #Oligodendrocyte Precursor Cells (OPC)
fplot_oli <- FeaturePlot(data.integrated, features = oli_gene, ncol=4)  #Oligodendrocytes
fplot_per <- FeaturePlot(data.integrated, features = per_gene, ncol=4)  #Pericytes
fplot_end <- FeaturePlot(data.integrated, features = end_gene, ncol=4)  #Endothelial Cells




#####figure 9 and 10 #####
oli_gene <- c( "MOG",	"MAG",	"PLP1",	 "CNP")
ast_gene <- c("AQP4",	 "GFAP",	 "FGFR3",	 "GJA1")
mic_gene <- c( "CSF1R",	 "P2RY12",	 "CD74",	"C1QA")
opc_gene <- c( "PDGFRA",	 "BCAN",	"CSPG4",	 "PTPRZ1")
inh_gene <- c("GAD1",	 "GAD2",	 "ADARB2",	 "LHX6")
per_gene <- c( "ATP1A2",	 "CLDN5",	 "DCN",	 "EPS8")
end_gene <- c("CLDN5",	 "ACTA2",	 "EPAS1",	 "FLT1")
exc_gene <- c("SLC17A7",	 "SATB2",	 "SNAP25",	 "FEZF2")

vplot_oli <- VlnPlot(data.integrated, features = oli_gene, pt.size = 0, assay = "RNA", ncol=4)
vplot_ast <- VlnPlot(data.integrated, features = ast_gene, pt.size = 0, assay = "RNA", ncol=4)  #Astrocytes 4
vplot_mic <- VlnPlot(data.integrated, features = mic_gene, pt.size = 0, assay = "RNA", ncol=4)  #Microglia 5
vplot_opc <- VlnPlot(data.integrated, features = opc_gene, pt.size = 0, assay = "RNA", ncol=4)  #Oligodendrocyte Precursor Cells (OPC) 7
vplot_inh <- VlnPlot(data.integrated, features = inh_gene, pt.size = 0, assay = "RNA", ncol=4)  #Inhibitory Neurons 3

vplot_oli/vplot_ast/vplot_mic/vplot_opc/vplot_inh+
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ")") &  
  theme(axis.title.x = element_blank())

vplot_per/vplot_end/vplot_exc+
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = ")") &  
  theme(axis.title.x = element_blank())



######annotation#####
new.cluster.ids <-c("Excitatory Neurons","Oligodendrocytes", "Inhibitory Neurons", "Astrocytes", "Microglia", "OPC")
names(new.cluster.ids) <- levels(data.integrated)
data.integrated <- RenameIdents(data.integrated, new.cluster.ids)
dimplot_labeled <- DimPlot(data.integrated, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()




#####stacked bar chart#####
region = c(rep("Amygdala",6), rep("Cerebellum",6), rep("Hippocampus",6), rep("Occipital", 6), rep("Parietal",6)) 
celltype= rep(c("Excitatory Neurons","Oligodendrocytes", "Inhibitory Neurons", "Astrocytes", "Microglia", "OPC"),5)
df = data.frame(region,celltype)
df['value'] <- NA

for (i in 1:6){
  df[i,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Amygdala"))
  df[i+6,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Cerebellum"))
  df[i+12,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Hippocampus"))
  df[i+18,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Occipital"))
  df[i+24,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Parietal"))
}


percent_region <- ggplot(df, aes(fill=celltype, y=value, x=region)) + 
  geom_bar(position="fill", stat="identity")

percent_celltype   <- ggplot(df, aes(fill=region, y=value, x=celltype)) + 
  geom_bar(position="fill", stat="identity")


#####pie chart #####
# columns = c("Amygdala", "Cerebellum", "Hippocampus", "Occipital", "Parietal") 
# rows= c("Oligodendrocytes", "Cluster 1", "Cluster 2", "Inhibitory Neurons", "Astrocytes", "Microglia", "Cluster 6", "OPC")
# 
# df = data.frame(matrix(nrow = length(rows), ncol = length(columns))) 
# colnames(df) = columns
# rownames(df) = rows
# 
# for (i in 1:8){
#   df[i,1] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Amygdala"))
#   df[i,2] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Cerebellum"))
#   df[i,3] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Hippocampus"))
#   df[i,4] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Occipital"))
#   df[i,5] = nrow(filter(data.integrated@meta.data, seurat_clusters == i-1, region == "Parietal"))
# }
#
# for (p in 1:8){
#   #work out %
#   pct <- round(100*df[p,]/sum(df[p,]))
#   # Draw pie chart
#   pie(as.numeric(df[p,]),
#       labels = paste(columns, sep = " ", pct, "%"), 
#       #col = rainbow(length(df[p,])), 
#       main = paste(rows[p], sep = " ",": % of each regions"))
# }



#######proportion of gene markers by region######

cell_num_exc = nrow(filter(data.integrated@meta.data, seurat_clusters == 0))
cell_num_oli = nrow(filter(data.integrated@meta.data, seurat_clusters == 1))
cell_num_inh = nrow(filter(data.integrated@meta.data, seurat_clusters == 2))
cell_num_ast = nrow(filter(data.integrated@meta.data, seurat_clusters == 3))
cell_num_mic = nrow(filter(data.integrated@meta.data, seurat_clusters == 4))
cell_num_opc = nrow(filter(data.integrated@meta.data, seurat_clusters == 5))

#example
PLP1_df <- as.data.frame(GetAssayData(data.integrated, slot = "data", assay = "integrated"))["PLP1",]
PLP1_df <- filter(as.data.frame(t(PLP1_df)),  PLP1 >=2)
PLP1_reads <- rownames(PLP1_df)

A <- length(grep("Amygdala", PLP1_reads, ignore.case=FALSE))
C <- length(grep("Cerebellum", PLP1_reads, ignore.case=FALSE))
H <- length(grep("Hippocampus", PLP1_reads, ignore.case=FALSE))
O <- length(grep("Occipital", PLP1_reads, ignore.case=FALSE))
P <- length(grep("Parietal", PLP1_reads, ignore.case=FALSE))




# save figures
# plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
# plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
# file.copy(from=plots.png.paths, to="results_dir")

#session.info()