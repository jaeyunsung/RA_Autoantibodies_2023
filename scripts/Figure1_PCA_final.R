library("dplyr")
library("ggplot2")
library("FactoMineR")
library("stats")
library("factoextra")
library("ggfortify")

quantile_df <- read.csv("../data/all_quantile_df.csv")

quantile_df$Sample_Status <- as.factor(quantile_df$Sample_Status)

#Why not just read these data? we saved them as a file.
negative_df <- filter(quantile_df, quantile_df$Sample_Status == 'ACPA_negative')
positive_df <- filter(quantile_df, quantile_df$Sample_Status == 'ACPA_positive')
control_df <- filter(quantile_df, quantile_df$Sample_Status == 'control')

col_nums = dim(positive_df)[2]
num_features1 <- ncol(quantile_df)-2
feature_names_only <- data.frame(colnames(positive_df[,1:num_features1]))

ANOVA_pvalues1 <- data.frame(matrix(data = 0, nrow = num_features1, ncol=2))
ANOVA_pvalues1[1:num_features1,1] <- feature_names_only
colnames(ANOVA_pvalues1) <- cbind("Autoantibody","Pvalue")

for (i in 1:num_features1) {
    one.way2 <- aov( lm(get(feature_names_only[i,1]) ~ get("Sample_Status"), quantile_df))
    #one.way2 <- aov( lm(get(feature_names_only[i,1]) ~ get("Sample_Status"), quantile_df), data = quantile_df
    p_val <- summary(one.way2)[[1]][1,5]
    ANOVA_pvalues1[i,2] <- p_val
}

ordered_anova_pvalues <- ANOVA_pvalues1[order(ANOVA_pvalues1$Pvalue),]

#LINE 92 - 141 can be reduced with this.
features_less_than_0.05 <- ordered_anova_pvalues[which(ordered_anova_pvalues["Pvalue"] < 0.05),]
num_sig_features <- dim(features_less_than_0.05)[1]

sig_feature_list <- features_less_than_0.05[,1]
sig_feature_df = quantile_df[,which(colnames(quantile_df) %in% sig_feature_list)]
sig_feature_df$Sample_Status <- quantile_df$Sample_Status
sig_feature_df$Case_ID <- quantile_df$Case_ID

#print(sig_feature_df)
sig_feature_table <- sig_feature_df[, 1:26]
#print(sig_feature_table)
pca2 <- prcomp(sig_feature_table, scale. = F)
PCA1 <- autoplot(pca2, data = sig_feature_df, colour = 'Sample_Status', 
         frame.type = 'norm')+ scale_color_manual(values = c("gray50","brown1", "steelblue3"))+ coord_fixed()
pdf("../output/PCA.pdf")
PCA1
dev.off()


