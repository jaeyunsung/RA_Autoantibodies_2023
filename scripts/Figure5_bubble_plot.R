
library("ggpubr")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("reshape2")

#read in the quantile normalized data
quantile_df <- read.csv("../data/quantile_for_CDAI_scatterplots.csv")
cdai_df2 <- read.csv("../data/patient_info_table_v3.csv", header = TRUE)



print(cdai_df2$Sample_ID == quantile_df$Sample_ID)


#get number of features
num_features1 <- ncol(quantile_df)-2
num_patients <- dim(quantile_df)[1]


#get feature names
feature_name_df <- colnames(quantile_df[1:num_features1])
#make a dataframe to hold the results
spearman_values_df <- as.data.frame(matrix(data = 0, nrow = num_features1, ncol = 3))
spearman_values_df[1:num_features1,1] <- feature_name_df
colnames(spearman_values_df) <- cbind("Autoantibody", "Spearman_value_all_RA", "P_value_all_RA")
#View(spearman_values_df)
#calculate spearman correlations

feature1 <- cdai_df2$CDAI
#View(feature1)
for (i in 1:num_features1) {
  feature2 <- quantile_df[1:num_patients, i]
  pearson_value <- cor.test(feature1, feature2, method = "spearman")[4]
  spearman_values_df[i,2] <- pearson_value
  pval <- cor.test(feature1, feature2, method = "spearman")[3]
  spearman_values_df[i,3] <- pval
}
#View(spearman_values_df)


#Now split ACPA+ and ACPA-  groups
ACPA_pos_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_positive_RA')
#View(ACPA_pos_df)
ACPA_neg_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_negative_RA')
#View(ACPA_neg_df)

ACPA_pos_CDAI <- filter(cdai_df2, cdai_df2$Study_group == 'ACPA_positive_RA')
#View(ACPA_pos_CDAI)
ACPA_neg_CDAI <- filter(cdai_df2, cdai_df2$Study_group == 'ACPA_negative_RA')
#View(ACPA_neg_CDAI)

ACPA_pos_df <- ACPA_pos_df[order(ACPA_pos_df$Sample_ID),]
#View(ACPA_pos_df)
ACPA_pos_CDAI <- ACPA_pos_CDAI[order(ACPA_pos_CDAI$Sample_ID),]
#View(ACPA_pos_CDAI)
print(ACPA_pos_df$Sample_ID == ACPA_pos_CDAI$Sample_ID)

#Now do the same for ACPA-negative
ACPA_neg_df <- ACPA_neg_df[order(ACPA_neg_df$Sample_ID),]
#View(ACPA_neg_df)
ACPA_neg_CDAI <- ACPA_neg_CDAI[order(ACPA_neg_CDAI$Sample_ID),]
#View(ACPA_neg_CDAI)
print(ACPA_neg_df$Sample_ID == ACPA_neg_CDAI$Sample_ID)

#get number of samples for ACPA+ and ACPA-
num_pos <- dim(ACPA_pos_CDAI)[1]
print("number of ACPA positive")
print(num_pos)
num_neg <- dim(ACPA_neg_CDAI)[1]
print("number of ACPA negative")
print(num_neg)


#Make a data frame to hold the values
ACPA_pos_spearman <- as.data.frame(matrix(data = 0, nrow = num_features1, ncol = 3))
ACPA_pos_spearman[1:num_features1,1] <- feature_name_df
colnames(ACPA_pos_spearman) <- cbind("Autoantibody", "Spearman_value_ACPA_pos", "P_value_ACPA_pos")
#View(ACPA_pos_spearman)
feature1 <- ACPA_pos_CDAI$CDAI
for (i in 1:num_features1) {
  feature2 <- ACPA_pos_df[1:num_pos, i]
  spearman_val <- cor.test(feature1, feature2, method = "spearman")[4]
  ACPA_pos_spearman[i,2] <- spearman_val
  p_val <- cor.test(feature1, feature2, method = "spearman")[3]
  ACPA_pos_spearman[i,3] <- p_val
}
#View(ACPA_pos_spearman)


#Now do the same for ACPA-negative
ACPA_neg_spearman <- as.data.frame(matrix(data = 0, nrow = num_features1, ncol = 3))
ACPA_neg_spearman[1:num_features1,1] <- feature_name_df
colnames(ACPA_neg_spearman) <- cbind("Autoantibody", "Spearman_value_ACPA_neg", "P_value_ACPA_neg")
#View(ACPA_neg_spearman)

feature1 <- ACPA_neg_CDAI$CDAI
for (i in 1:num_features1) {
  feature2 <- ACPA_neg_df[1:num_neg, i]
  spearman_val <- cor.test(feature1, feature2, method = "spearman")[4]
  ACPA_neg_spearman[i,2] <- spearman_val
  p_val <- cor.test(feature1, feature2, method = "spearman")[3]
  ACPA_neg_spearman[i,3] <- p_val
}
#View(ACPA_neg_spearman)


#combine the results of the three datasets
all_spearman_results <- cbind(spearman_values_df, ACPA_pos_spearman[1:num_features1,2:3], ACPA_neg_spearman[1:num_features1,2:3])
#View(all_spearman_results)

#now get those that have a p value less than 0.01
all_spearman_sig_pvals <- filter(all_spearman_results, all_spearman_results$`P_value_all_RA` < 0.01 | all_spearman_results$`P_value_ACPA_pos` < 0.01 | all_spearman_results$`P_value_ACPA_neg` < 0.01)
#View(all_spearman_sig_pvals)


#now remove those that have a spearman value less than 0.4 and greater than -0.4
all_significant_spearman_and_pvals <- filter(all_spearman_sig_pvals, (all_spearman_sig_pvals$Spearman_value_all_RA > 0.4 & all_spearman_sig_pvals$`P_value_all_RA` < 0.01) | (all_spearman_sig_pvals$Spearman_value_all_RA < -0.4 & all_spearman_sig_pvals$`P_value_all_RA` < 0.01)| (all_spearman_sig_pvals$Spearman_value_ACPA_pos > 0.4 & all_spearman_sig_pvals$`P_value_ACPA_pos` < 0.01 ) | 
(all_spearman_sig_pvals$Spearman_value_ACPA_pos < -0.4 & all_spearman_sig_pvals$`P_value_ACPA_pos` < 0.01) | 
(all_spearman_sig_pvals$Spearman_value_ACPA_neg > 0.4 & all_spearman_sig_pvals$`P_value_ACPA_neg` < 0.01) | (all_spearman_sig_pvals$Spearman_value_ACPA_neg < -0.4 & all_spearman_sig_pvals$`P_value_ACPA_neg` < 0.01))
#View(all_significant_spearman_and_pvals)


temp_df2 <- data.frame(all_significant_spearman_and_pvals)
temp_df2$log_all_RA <- -log10(temp_df2$P_value_all_RA)
temp_df2$log_ACPA_pos <- -log10(temp_df2$P_value_ACPA_pos)
temp_df2$log_ACPA_neg <- -log10(temp_df2$P_value_ACPA_neg)
temp_df2$ACPA_neg <- 1.5
temp_df2$ACPA_pos <- 2
temp_df2$All_RA <- 2.5

#Order of samples: ACPA+, ACPA-  , ALL RA
bubble_plot1 <- ggplot()+geom_point(data = temp_df2, aes(y = Autoantibody, x=ACPA_pos, colour = Spearman_value_ACPA_pos, size=log_ACPA_pos )) + geom_point(data = temp_df2, aes(y = Autoantibody, x=ACPA_neg, colour = Spearman_value_ACPA_neg, size=log_ACPA_neg ))+geom_point(data = temp_df2, aes(y = Autoantibody, x=All_RA, colour = Spearman_value_all_RA, size=log_all_RA ))+ xlim(1,3) + 
  scale_colour_gradient2(low = "steelblue", mid = "white",high = "red", breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5))+scale_size_continuous(breaks = c((-log10(0.001)),(-log10(0.01)),(-log10(0.1)), (-log10(0.2))))+ theme(axis.text.x=element_blank(),axis.text = element_text(face="bold"))+xlab("RA Patients")
bubble_plot1

pdf('../output/Bubble_plot.pdf')
bubble_plot1
dev.off()








