## -----------------------------------------------------------------------------
library("ggpubr")
library("ggplot2")


## -----------------------------------------------------------------------------
#cdai_df2 <- read.csv("/Users/Kevin/Desktop/Grad school docs and past classes/New_RA_Project_code/Github_code/patient_info_table_v3.csv", header = TRUE)

cdai_df2 <- read.csv("../data/patient_info_table_v3.csv", header = TRUE)

#View(cdai_df2)


## -----------------------------------------------------------------------------
quantile_df <- read.csv("../data/quantile_for_CDAI_scatterplots.csv")
#View(quantile_df)


## -----------------------------------------------------------------------------
dim(quantile_df)
dim(cdai_df2)


## -----------------------------------------------------------------------------
cdai_df2$slideID == quantile_df$Case_ID


## -----------------------------------------------------------------------------
num_features1 <- ncol(quantile_df)-2
#num_features1


## -----------------------------------------------------------------------------
feature_name_df <- colnames(quantile_df[1:num_features1])
#View(feature_name_df)


## -----------------------------------------------------------------------------
spearman_values_df <- as.data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
spearman_values_df[1:1622,1] <- feature_name_df
colnames(spearman_values_df) <- cbind("Autoantibody", "Spearman_value")
#View(spearman_values_df)


## -----------------------------------------------------------------------------
feature1 <- cdai_df2$cdai
#View(feature1)
for (i in 1:num_features1) {
  feature2 <- quantile_df[1:56, i]
  pearson_value <- cor.test(feature1, feature2, method = "spearman")[4]
  spearman_values_df[i,2] <- pearson_value
}
#View(spearman_values_df)



## -----------------------------------------------------------------------------
ordered_spearman_values_df <- spearman_values_df[order(-spearman_values_df$Spearman_value),]
#View(ordered_spearman_values_df)


## -----------------------------------------------------------------------------
features_to_plot_spearman <- data.frame(matrix(data = 0, nrow = 56, ncol = 3))
colnames(features_to_plot_spearman) <- cbind("CDAI","CISH", "Sample_Status")
features_to_plot_spearman$CDAI <- cdai_df2$cdai
features_to_plot_spearman$CISH <- quantile_df$CISH
features_to_plot_spearman$Sample_Status <- quantile_df$Sample_Status



## -----------------------------------------------------------------------------
scatter_plot <- ggplot(data = features_to_plot_spearman, aes(x = CDAI, y = CISH))+ geom_point(aes(color=Sample_Status))+ geom_smooth(method="lm", se=FALSE)+stat_cor(method = "spearman")+ scale_colour_manual(values = c("gray50","brown1"))+ coord_fixed()
#


## -----------------------------------------------------------------------------

pdf('../output/Spearman_plot.pdf')
scatter_plot
dev.off()

## -----------------------------------------------------------------------------

