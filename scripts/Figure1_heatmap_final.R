
 
library(pheatmap)

 

#Building a heatmap of all features for all 92 samples.
#Read in all samples
 
all_quantile <- read.csv("../data/all_quantile_df.csv")
#View(all_quantile_df)
 

#Check dimensions
 
dim(all_quantile)
 

#get number of features and samples
 
num_features <- ncol(all_quantile)-2
#num_features
num_samples <- nrow(all_quantile)
#num_samples
 

#Remove the column of Case IDs and sample status. Make the row names the Case ID
 
all_quantile_df <- data.frame(matrix(data = 0, nrow=num_samples, ncol=num_features))
rownames(all_quantile_df) <- all_quantile$Case_ID
colnames(all_quantile_df) <- colnames(all_quantile[1,1:num_features])
all_quantile_df[1:num_samples,1:num_features] <- all_quantile[1:num_samples,1:num_features]
#View(all_quantile_df)

 

#Make a vector of the sample status for the heatmap
 
colors_df <- data.frame(matrix(data = 0, nrow = 92, ncol=1))
colnames(colors_df) <- c("Sample_Status")
colors_df[1:92,1] <- all_quantile$Sample_Status
rownames(colors_df) <- all_quantile$Case_ID
#View(colors_df)
 

#Define the breaks and the colors
 
library(RColorBrewer)
break_list = seq(-1.0,1.0, by=0.05)
color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(break_list))
 

#Now use z score
 
quantile_zscore_df <- scale(all_quantile_df, center = TRUE, scale = TRUE)
#View(quantile_zscore_df)
quantile_zscore_transposed <- t(quantile_zscore_df)
#View(quantile_zscore_transposed)
 

#Plot
 
heatmap1 <- pheatmap(quantile_zscore_transposed, annotation = colors_df, show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = TRUE,  breaks = break_list, color = color)
pdf("../output/heatmap.pdf")
heatmap1
dev.off()





