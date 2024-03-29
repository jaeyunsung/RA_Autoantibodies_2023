library("dplyr")
library("ggplot2")
library("FactoMineR")
library("stats")
library("factoextra")
library("ggfortify")

quantile_df <- read.csv("../data/all_quantile_df.csv")

#quantile_df$Study_group <- as.factor(quantile_df$Study_group)


ACPA_neg_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_negative_RA')
ACPA_pos_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_positive_RA')
control_df <- filter(quantile_df, quantile_df$Study_group == 'control')

quantile_df <- rbind(ACPA_neg_df,ACPA_pos_df,control_df)
#print(quantile_df$Study_group)
#print(dim(quantile_df))

col_nums = dim(ACPA_pos_df)[2]
num_features1 <- ncol(quantile_df)-2
feature_names_only <- data.frame(colnames(ACPA_pos_df[,1:num_features1]))

ANOVA_pvalues1 <- data.frame(matrix(data = 0, nrow = num_features1, ncol=2))
ANOVA_pvalues1[1:num_features1,1] <- feature_names_only
colnames(ANOVA_pvalues1) <- cbind("Autoantibody","Pvalue")

for (i in 1:num_features1) {
    one.way2 <- aov( lm(get(feature_names_only[i,1]) ~ get("Study_group"), quantile_df))
    p_val <- summary(one.way2)[[1]][1,5]
    ANOVA_pvalues1[i,2] <- p_val
}

ordered_anova_pvalues <- ANOVA_pvalues1[order(ANOVA_pvalues1$Pvalue),]


features_less_than_0.05 <- ordered_anova_pvalues[which(ordered_anova_pvalues["Pvalue"] < 0.05),]
num_sig_features <- dim(features_less_than_0.05)[1]

sig_feature_list <- features_less_than_0.05[,1]
sig_feature_df = quantile_df[,which(colnames(quantile_df) %in% sig_feature_list)]
sig_feature_df$Study_group <- quantile_df$Study_group
sig_feature_df$Sample_ID <- quantile_df$Sample_ID

#print(sig_feature_df)
sig_feature_table <- sig_feature_df[, 1:26]
#print(sig_feature_table)
pca2 <- prcomp(sig_feature_table, scale. = F)
PCA1 <- autoplot(pca2, data = sig_feature_df, colour = 'Study_group', 
         frame.type = 'norm')+ scale_color_manual(values = c("gray50","brown1", "steelblue3"))+ coord_fixed()
pdf("../output/PCA.pdf")
PCA1
dev.off()


