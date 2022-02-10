#First I need to remove
#outlier 1: Case028292 ACPA Negative
#outlier 2: Case028399 ACPA Positive
#for b cell therapy

#Second, remove the 2 outliers from seronegative
#before doing quantile normalization
#Case028344
#Case028311

#Third, remove any that have a CDAI score of 0, they are NAs
#Case028327 ACPA_negative
#Case028276 ACPA_negative
#Case028368 ACPA_negative
#Case028369 ACPA_positive
#Case028281 ACPA_positive
#Case028386 ACPA_positive

#Now there are 66 - 10 = 56 samples with CDAI scores to use
#Case028399

 
library("ggpubr")
library("ggplot2")
 

#Read in Clinical data
 
cdai_df2 <- read.csv("../RA_Autoantibodies_2022/patient_info_table_v3.csv", header = TRUE)
#View(cdai_df2)
 

#Read in autoantibody data
 
quantile_df <- read.csv("../RA_Autoantibodies_2022/quantile_for_CDAI_scatterplots.csv")
#View(quantile_df)
 

#Check how many samples are in each dataset
 
dim(quantile_df)
dim(cdai_df2)
 

#Make sure that the case ids are in the same order
 
#cdai_df2$slideID == quantile_df$Case_ID
 

#Get number of features
 
num_features1 <- ncol(quantile_df)-2
#num_features1
 

#Get feature names
 
feature_name_df <- colnames(quantile_df[1:num_features1])
#View(feature_name_df)
 

#Make a dataframe for the results
 
spearman_values_df <- as.data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
spearman_values_df[1:1622,1] <- feature_name_df
colnames(spearman_values_df) <- cbind("Autoantibody", "Spearman_value")
#View(spearman_values_df)
 

#Now do the spearman correlations
 
feature1 <- cdai_df2$cdai
#View(feature1)
for (i in 1:num_features1) {
  feature2 <- quantile_df[1:56, i]
  pearson_value <- cor.test(feature1, feature2, method = "spearman")[4]
  spearman_values_df[i,2] <- pearson_value
}
#View(spearman_values_df)

 

#Order them and plot the autoantibody with the highest Spearman value
 
ordered_spearman_values_df <- spearman_values_df[order(-spearman_values_df$Spearman_value),]
#View(ordered_spearman_values_df)
 

#CISH has the largest spearman
 
features_to_plot_spearman <- data.frame(matrix(data = 0, nrow = 56, ncol = 3))
colnames(features_to_plot_spearman) <- cbind("CDAI","CISH", "Sample_Status")
features_to_plot_spearman$CDAI <- cdai_df2$cdai
features_to_plot_spearman$CISH <- quantile_df$CISH
features_to_plot_spearman$Sample_Status <- quantile_df$Sample_Status
#View(features_to_plot_spearman)
 

#Add Sample status to this plot
 
scatter_plot <- ggplot(data = features_to_plot_spearman, aes(x = CDAI, y = CISH))+ geom_point(aes(color=Sample_Status))+ geom_smooth(method="lm", se=FALSE)+stat_cor(method = "spearman")+ scale_colour_manual(values = c("gray50","brown1"))+ coord_fixed()

pdf("../RA_Autoantibodies_2022/output/Spearman_plot.pdf")
scatter_plot
dev.off() 