# Load Libraries ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(data.table)
library(readxl)
library(viridis)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(patchwork)
# Import and Filter -------------------------------------------------------


#10x Data importation from cell ranger output files, takes barcodes.tsv, features.tsc, matrix.mtx.
hpf36_OMP_RFP <- Read10X(data.dir = "36hrMergedCount")
hpf48_OMP_RFP <- Read10X(data.dir = "48hrMergedCount")
hpf72_OMP_RFP <- Read10X(data.dir = "72hrMergedCount")

# Initialize the Seurat object with the raw (non-normalized data).
hpf36_OMP_RFP <- CreateSeuratObject(counts = hpf36_OMP_RFP, project = "36hpf_OMP_RFP", min.cells = 2, min.features = 100)
hpf48_OMP_RFP <- CreateSeuratObject(counts = hpf48_OMP_RFP, project = "48hpf_OMP_RFP", min.cells = 2, min.features = 100)
hpf72_OMP_RFP <- CreateSeuratObject(counts = hpf72_OMP_RFP, project = "72hpf_OMP_RFP", min.cells = 2, min.features = 100)

combi_OMP_RFP <- merge(hpf36_OMP_RFP, y = c(hpf48_OMP_RFP, hpf72_OMP_RFP), add.cell.ids = c("36hpf", "48hpf", "72hpf"))
combi_OMP_RFP
combi_OMP_RFP@meta.data
levels(combi_OMP_RFP)

#Start quality control analysis. First quantify mitochondrial expression. 
combi_OMP_RFP[["percent.mt"]] <- PercentageFeatureSet(combi_OMP_RFP, pattern = "^mt-")
combi_OMP_RFP[["percent.OR"]] <- PercentageFeatureSet(combi_OMP_RFP, pattern = "^or1")

#visualize qc data
VlnPlot(combi_OMP_RFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "percent.OR")
CombinePlots(plots = list(plot1, plot2))
remove(plot1,plot2)
#Stats before filtering
combi_OMP_RFP
#Filter
combi_OMP_RFP <- subset(combi_OMP_RFP, subset = nFeature_RNA > 800 & nFeature_RNA < 7000 & percent.mt < 10 & nCount_RNA < 50000)
#stats after filtering
combi_OMP_RFP
plot1 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
remove(plot1,plot2)
VlnPlot(combi_OMP_RFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "percent.OR")
VlnPlot(combi_OMP_RFP, features = c("nFeature_RNA", "nCount_RNA", "percent.OR"), ncol = 3)
combi_OMP_RFP@meta.data

# Normalization -----------------------------------------------------------
#Using SCTransform to normalize the data
combi_OMP_RFP<- SCTransform(combi_OMP_RFP,return.only.var.genes = FALSE, verbose = TRUE)

#features post normalization
VlnPlot(combi_OMP_RFP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combi_OMP_RFP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
remove(plot1,plot2)
#Export SCT normalized Data to Excel. Uses Data.table package. Add date to file name if wanted. 
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP[["SCT"]]@counts))
fwrite(x = data_to_write_out, row.names = TRUE, file = "SCTransAll_combi_Outfile_10012021.csv")
remove(data_to_write_out)
saveRDS(combined_OSNs, file = "combined_OSNs_NormalizedSCT.rds")

# OR naming ---------------------------------------------------------------
#Imports excel files from previously written out data that has been edited.
#I manually removed all genes that were not ORs, I wanted this excel file to 
#be able to manually check some parameters outside of R
#Read in data
OR_df <- readxl::read_xlsx("SCTransAll_combi_Outfile_OR_10012021.xlsx")

#Lets deal with the OR ratios and identities first
#Convert to data.table, set the column of genes as rowname
OR_dt <- as.data.table(OR_df[-1])
rownames(OR_dt) <- OR_df[[1]]

# filter out cells that have no OR genes
OR_all_zero_cols <- which(sapply(OR_dt, function(x) sum(x) == 0))
OR_dt[, (OR_all_zero_cols) := NULL]

# compute information about prominent genes
OR_top2s <- lapply(OR_dt, function(x) {
  names(x) <- rownames(OR_dt)
  top2 <- tail(sort(x), 2)
  list(first = names(top2)[2], second = names(top2)[1], ratio = top2[1]/top2[2])
})

gene_prominence_data <- rbindlist(OR_top2s)[, gene := names(OR_top2s)][ratio <= 1]

# pull cells with prominent genes and the name of the prominent gene
prominent_genes <- setNames(gene_prominence_data$first, gene_prominence_data$gene)
prominent_genes_ratios <- setNames(gene_prominence_data$ratio, gene_prominence_data$gene)

# merge with other cells that don't have a prominent gene (default value NA)
prominent_genes_metadata <- setNames(rep(NA_character_, ncol(OR_df) - 1), colnames(OR_df)[-1])
prominent_genes_metadata[names(prominent_genes)] <- prominent_genes

#attempting ratio
prominent_genes_ratios_metadata <- setNames(rep(NA_character_, ncol(OR_df) - 1), colnames(OR_df)[-1])
prominent_genes_ratios_metadata[names(prominent_genes_ratios)] <- prominent_genes_ratios

#Copy 0 OR ratio (only cells with 1 OR expressed)
combi_OMP_RFP <- AddMetaData(
  object = combi_OMP_RFP,
  metadata = prominent_genes_metadata,
  col.name = 'prominent.OR'
)
combi_OMP_RFP <- AddMetaData(
  object = combi_OMP_RFP,
  metadata = as.numeric(prominent_genes_ratios_metadata),
  col.name = 'prominent.OR.ratio'
)
# check if it's been added
head(combi_OMP_RFP@meta.data)
# should be accessible w/ double brackets
combi_OMP_RFP [["prominent.OR"]]
combi_OMP_RFP [["prominent.OR.ratio"]]

#Write out data
#export Meta data showing prominent OR per cell and OR ratio
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP@meta.data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "10122021_combi_OMP_RFP_metadata.csv")
remove(data_to_write_out)
#export co-expression OR information showing top 2 ORs and ratio
data_to_write_out <- as.data.frame(as.matrix(gene_prominence_data))
fwrite(x = data_to_write_out, row.names = TRUE, file = "10122021_combi_OMP_RFP_OR_prominence_data.csv")
remove(data_to_write_out)

#cleanup
remove(gene_prominence_data, OR_df, OR_dt, OR_top2s, OR_all_zero_cols, prominent_genes, prominent_genes_metadata, prominent_genes_ratios, prominent_genes_ratios_metadata)

# OR Subsets --------------------------------------------------------------
#Create alternate object for OMP_RFP based on assigned ORs, only cells with called ORs will remain

#OR_combi_OMP_RFP <- subset(combi_OMP_RFP, subset = prominent.OR != "N/A")

#OR called if 5x expression of second OR
#pt2_combi_OMP_RFP <- subset(combi_OMP_RFP, subset = prominent.OR != "N/A" & prominent.OR.ratio <= 0.2)

# Processing Code ---------------------------------------------------------
#This is the standard data processing workflow I settled on using 
#It can be used for any data subset you want to analyse

#Process all all OSNs
remove(tmp,all.genes, top30)
tmp <- combi_OMP_RFP
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
top30 <- head(VariableFeatures(tmp), 30)
plot1 <- VariableFeaturePlot(tmp)
LabelPoints(plot = plot1, points = top30, repel = TRUE)
remove(plot1)
#Linear transformation for data scaling
all.genes <- rownames(tmp)
tmp <- ScaleData(tmp, features = all.genes)
tmp<- RunPCA(tmp, features = VariableFeatures(object = tmp))
print(tmp[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(tmp, dims = 1:2, reduction = "pca")
DimPlot(tmp, reduction = "pca")
DimHeatmap(tmp, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(tmp, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(tmp)
#Clustering
tmp <- FindNeighbors(tmp, dims = 1:10)
tmp <- FindClusters(tmp, resolution = 0.5)
###head(Idents(tmp, 5))
#UMAP Dimensionality Reduction
tmp <- RunUMAP(tmp, dims = 1:10)
DimPlot(tmp, reduction = "umap", label = TRUE, label.color = 'white') +DarkTheme()
DimPlot(tmp, reduction = "umap", label = TRUE)

# Output ------------------------------------------------------------------
#reassign temp to variable
combi_OMP_RFP<-tmp
remove(tmp)
levels(combi_OMP_RFP)

#Save RDS backup and Workspace if wanted
saveRDS(combi_OMP_RFP, file = "combined_OSNs_DataProcessed.rds")
save.image()

# Data Visualization Figures -----------------------------------------------

##Figure 1
# Visualize UMAPs and Clustering
remove(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
p0<-FeaturePlot(combi_OMP_RFP, features = "ompb") * scale_colour_viridis()* DarkTheme() + ggtitle('Placeholder for OMP :RFP Image', )
p0
p1<-DimPlot(combi_OMP_RFP, reduction = "umap", label = TRUE, pt.size = .1,label.color = "White",  label.box = TRUE, repel = TRUE)  + theme(plot.title = element_text(hjust = 0.5)) * DarkTheme()
p1
p2<-DimPlot(combi_OMP_RFP, reduction = "umap", label = TRUE, split.by = 'orig.ident', label.size = 0, label.color = 'white', pt.size = .1, shuffle = TRUE) + DarkTheme() + NoLegend()
p2

# Define the layout design 
combined_plot <- (p0 | p1) /
  p2 + 
  plot_layout(heights = c(1, 1))
combined_plot

#Redefining Clusters by maturity
Combined_clusters <- combi_OMP_RFP

# Change seurat_clusters to the cluster column name
new_cluster_ids <- Combined_clusters@meta.data$seurat_clusters 
new_cluster_ids <- as.character(Combined_clusters@meta.data$seurat_clusters)

# Assign new group names
new_cluster_ids[new_cluster_ids %in% c(3,5,8,4,10,11,13,14)] <- 'Immature'
new_cluster_ids[new_cluster_ids %in% c(2,0)] <- 'Early Mature'
new_cluster_ids[!new_cluster_ids %in% c('Immature', 'Early Mature',7)] <- 'Mature'

# Update the metadata with the new clusters
Combined_clusters@meta.data$new_clusters <- new_cluster_ids

# Update the identity classes
Combined_clusters <- SetIdent(Combined_clusters, value = 'new_clusters')

#Visualize new clusters IDS Including cluster 7
DimPlot(Combined_clusters, reduction = "umap", label = TRUE, pt.size = .1, label.color = "White",  label.box = TRUE, repel = TRUE)  + theme(plot.title = element_text(hjust = 0.5)) * DarkTheme()

#Remake without cluster 7
cells_to_keep <- WhichCells(Combined_clusters, expression = seurat_clusters != 7)

# Subset the Seurat object to keep only the cells you want
Combined_clusters <- subset(Combined_clusters, cells = cells_to_keep)

# seurat_obj_subset now contains your data minus the cells from cluster 7
DimPlot(Combined_clusters, reduction = "umap", label = TRUE, pt.size = .1,label.color = "White",  label.box = TRUE, repel = TRUE)  + theme(plot.title = element_text(hjust = 0.5)) * DarkTheme()

# Set the factor levels in the order you want
Combined_clusters$new_clusters <- factor(Combined_clusters$new_clusters, levels = c("Immature", "Early Mature", "Mature"))

maturity_markers <- c("neurod1", "ascl1a", "gap43","gng8", "ompb", "vamp2","elf3", "malb", "krt4","cfl1l", "cldne", "cfap157")

#Dotplot visualizing maturity marker by maturity groups
DotPlot(Combined_clusters, features = maturity_markers, group.by = "new_clusters", dot.scale = 20) + scale_color_viridis() + DarkTheme()

##Supplement 1, Maturity markers on all cluster including 7
FeaturePlot(combi_OMP_RFP, features = maturity_markers, ncol = 3) * scale_colour_viridis()* DarkTheme()

##Supplement 2 UMAPS of Immune Genes
Cluster7_Immune_genes <- c("apoc1", "tln1", "cd83","ccl35.1", "cxcr3.3", "tnfa", "tnfb","ptprc", "marco")
FeaturePlot(combi_OMP_RFP, features = Cluster7_Immune_genes, ncol = 3) * scale_colour_viridis()* DarkTheme()


# DEG Testing -------------------------------------------------------------
#Deseq2 Differentially Expressed Genes by cluster
# Ensure you're using the RNA assay for Deseq2
DefaultAssay(combi_OMP_RFP) <- "RNA"  

combi_OMP_RFP_15_Deseq2<- FindMarkers(combi_OMP_RFP, ident.1 = 15, test.use = "DESeq2")
combi_OMP_RFP_15_Deseq2
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP_9_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "combi_OMP_RFP_15_Deseq2_10072021.csv") 

combi_OMP_RFP_9v15_Deseq2<- FindMarkers(combi_OMP_RFP, ident.1 = 9, ident.2 = 15, test.use = "DESeq2")
combi_OMP_RFP_9v15_Deseq2
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP_9v15_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "combi_OMP_RFP_9v15_Deseq2_10072021.csv") 

combi_OMP_RFP_9v1_Deseq2<- FindMarkers(combi_OMP_RFP, ident.1 = 9, ident.2 = 1, test.use = "DESeq2")
combi_OMP_RFP_9v1_Deseq2
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP_9v1_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "combi_OMP_RFP_9v1_Deseq2_10072021.csv") 

combi_OMP_RFP_0a2v5a1_Deseq2<- FindMarkers(combi_OMP_RFP, ident.1 = c(0,2), ident.2 = c(5,1), test.use = "DESeq2")
combi_OMP_RFP_0a2v5a1_Deseq2
data_to_write_out <- as.data.frame(as.matrix(combi_OMP_RFP_0a2v5a1_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "combi_OMP_RFP_0a2v5a1_Deseq2_10072021.csv") 

#Counting ORs Export to excel
combi_pt2_ORs <- pt2_combi_OMP_RFP$prominent.OR
data_to_write_out <- as.data.frame(as.matrix(combi_pt2_ORs))
fwrite(x = data_to_write_out, row.names = TRUE, file = "combi_pt2_ORs_10132021.csv")

#To perform DESEQ2 analysis by OR Clade
#Define OR clades
CladeA_ORs <- c("or101-1", "or102-3", "or102-4" , "or102-5", "or103-1", "or103-2","or103-5", "or104-1", "or104-2", 
                "or106-1", "or106-2" , "or106-3", "or106-4", "or106-6", "or106-7", "or106-8", "or106-9", "or106-10", "or106-11", 
                "or107-1", "or108-1", "or108-2", "or108-3", "or109-1", "or109-2", "or109-7", "or109-11", "or109-13",
                "or110-2", "or111-1", "or111-2", "or111-3", "or111-4", "or111-5", "or111-6", "or111-7", "or111-8", "or111-9", "or111-10", "or111-11", "or112-1")

CladeB_ORs <- c("or115-1", "or115-2", "or115-5", "or115-6", "or115-7", "or115-9", "or115-10", "or115-11", "or115-12", "or115-15",
                "or116-1", "or116-2", "or117-1", "or118-2", "or119-2", "or120-1", "or121-1", "or122-1", "or122-2", "or123-1", "or124-1", 
                "or125-1", "or125-2", "or125-3", "or125-4", "or125-5", "or125-6", "or125-7", "or125-8", "or126-1", "or126-2", "or126-3", "or126-4", "or126-5",
                "or127-1", "or128-1", "or128-2", "or128-3", "or128-4", "or128-5", "or128-6", "or128-8", "or128-9", "or128-10")

CladeC_ORs <- c("or129-1", "or130-1", "or130-1.1", "or131-1", "or131-2", "or132-1", "or132-2", "or132-2.1", "or132-4", "or132-5", 
                "or133-1", "or133-2", "or133-4", "or133-5", "or133-6", "or133-7", "or133-9", "or134-1", "or135-1", "or136-1",
                "or137-2", "or137-3", "or137-4", "or137-6", "or137-7", "or137-8", "or137-9")

#subset OSNs which have a confidently assigned OR 
#remove cluster 7 from the main data set
cells_to_keep <- WhichCells(combi_OMP_RFP, expression = seurat_clusters != 7)

# Subset the Seurat object to keep only the cells you want
combi_OMP_RFP_no7 <- subset(combi_OMP_RFP, cells = cells_to_keep)

#Only keep OSNs in which an OR can confidently be assigned
pt2_combi_OMP_RFP <- subset(combi_OMP_RFP_no7, subset = prominent.OR != "N/A" & prominent.OR.ratio <= 0.2)

#Subset OR_OMP_RFP by clades, set new idents per subset, merge together, check output. 
CladeA_pt2_combi_OMP_RFP <- subset(x = pt2_combi_OMP_RFP, subset = prominent.OR %in% CladeA_ORs)
CladeB_pt2_combi_OMP_RFP <- subset(x = pt2_combi_OMP_RFP, subset = prominent.OR %in% CladeB_ORs)
CladeC_pt2_combi_OMP_RFP <- subset(x = pt2_combi_OMP_RFP, subset = prominent.OR %in% CladeC_ORs)

#Add meta data for each OR clade and then set idents for the cells in that clade
CladeA_pt2_combi_OMP_RFP<- AddMetaData(object = CladeA_pt2_combi_OMP_RFP, metadata = "CladeA", col.name = 'OR.Clade')
Idents(CladeA_pt2_combi_OMP_RFP) <- CladeA_pt2_combi_OMP_RFP$OR.Clade
CladeB_pt2_combi_OMP_RFP<- AddMetaData(object = CladeB_pt2_combi_OMP_RFP, metadata = "CladeB", col.name = 'OR.Clade')
Idents(CladeB_pt2_combi_OMP_RFP) <- CladeB_pt2_combi_OMP_RFP$OR.Clade
CladeC_pt2_combi_OMP_RFP<- AddMetaData(object = CladeC_pt2_combi_OMP_RFP, metadata = "CladeC", col.name = 'OR.Clade')
Idents(CladeC_pt2_combi_OMP_RFP) <- CladeC_pt2_combi_OMP_RFP$OR.Clade

#Merge the 3 subsets back together with there new idents

CladesABC_combi_pt2_Merged <- merge(x= CladeA_pt2_combi_OMP_RFP, y = list (CladeB_pt2_combi_OMP_RFP, CladeC_pt2_combi_OMP_RFP))

#Check levels of merged object (Should have three levels A,B and C)
levels(CladesABC_combi_pt2_Merged)

#You can either prep the SCT for find markers or remove the SCT assay.
# I remove it since it is incompatible with DESeq2
#PrepSCTFindMarkers(CladesABC_combi_pt2_Merged, assay = "SCT", verbose = TRUE)

#Deseq2 comparing OR Clades
DefaultAssay(CladesABC_combi_pt2_Merged) <- "RNA"
testCladesABC_combi_pt2_Merged <- CladesABC_combi_pt2_Merged[["SCT"]] <- NULL

#A vs B
CladesAvsB_combi_pt2_Merged_Deseq2<- FindMarkers(CladesABC_combi_pt2_Merged, ident.1 = "CladeA", ident.2 = "CladeB", test.use = "DESeq2")
CladesAvsB_combi_pt2_Merged_Deseq2
data_to_write_out <- as.data.frame(as.matrix(CladesAvsB_combi_pt2_Merged_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "CladesAvsB_combi_pt2_10132021.csv")
remove(data_to_write_out)
#AvsC
CladesAvsC_combi_pt2_Merged_Deseq2<- FindMarkers(CladesABC_combi_pt2_Merged, ident.1 = "CladeA", ident.2 = "CladeC", test.use = "DESeq2")
CladesAvsC_combi_pt2_Merged_Deseq2
data_to_write_out <- as.data.frame(as.matrix(CladesAvsC_combi_pt2_Merged_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "CladesAvsC_combi_pt2_10132021.csv")
remove(data_to_write_out)
#BvsC
CladesBvsC_combi_pt2_Merged_Deseq2<- FindMarkers(CladesABC_combi_pt2_Merged, ident.1 = "CladeB", ident.2 = "CladeC", test.use = "DESeq2")
CladesBvsC_combi_pt2_Merged_Deseq2
data_to_write_out <- as.data.frame(as.matrix(CladesBvsC_combi_pt2_Merged_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "CladesBvsC_combi_pt2_10132021.csv")
remove(data_to_write_out)

# AB vs C
CladesABvsC_combi_pt2_Merged_Deseq2<- FindMarkers(CladesABC_combi_pt2_Merged, ident.1 = c("CladeA","CladeB"), ident.2 = "CladeC", test.use = "DESeq2")
CladesABvsC_combi_pt2_Merged_Deseq2
data_to_write_out <- as.data.frame(as.matrix(CladesABvsC_combi_pt2_Merged_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "CladesABvsC_combi_pt2_10132021.csv")
remove(data_to_write_out)

# B vs AC
CladesBvsAC_combi_pt2_Merged_Deseq2<- FindMarkers(CladesABC_combi_pt2_Merged, ident.1 = c("CladeB"), ident.2 = c("CladeA", "CladeC"), test.use = "DESeq2")
CladesBvsAC_combi_pt2_Merged_Deseq2
data_to_write_out <- as.data.frame(as.matrix(CladesBvsAC_combi_pt2_Merged_Deseq2))
fwrite(x = data_to_write_out, row.names = TRUE, file = "CladesBvsAC_combi_pt2_10132021.csv")
remove(data_to_write_out)