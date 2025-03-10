#Load Libraries
library(pheatmap)
#OSNs expressing multiple ORs were exported as csv files for further analysis. 
#OSNs can be read in seperately and analyzed. 
#Read CSV File (Rows: genes, Columns: samples)
geneCounts <- read.csv("combined_OSNs_All_OSNs_with_multiple_ORS_ORs_only_ORDERED.csv", header = TRUE, row.names = 1)

##Calculate a correlation matrix. Peason is deafault method I used spearman. 
cor.matrix <- cor(t(geneCounts), method = "spearman")
##Visualization (heatmap)
pheatmap(cor.matrix, cluster_rows = FALSE, legend=TRUE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col =  0, fontsize_row = 8, fontsize_col = 8, breaks = seq(0, 1, length.out = 100))

#This was used to generat heatmaps for individual subfamilies.  
SFgeneCounts <- read.csv("combined_OSNs_All_OSNs_with_multiple_ORS_ORs_only_111sf.csv", header = TRUE, row.names = 1)
cor.matrix <- cor(t(SFgeneCounts), method = "spearman")
# Visualization with legend removed and font size increased.
pheatmap(cor.matrix, cluster_rows = FALSE, legend=FALSE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col =  0, fontsize_row = 12, fontsize_col = 12, breaks = seq(0, 1, length.out = 100))












