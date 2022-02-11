
 
library("dplyr")
library("ggplot2")
 

#Make a function to get the means from a dataset for all features
 
feature_means <- function(results_df,data_df,num_features){
  num_row <- nrow(data_df)
  feature1 <- matrix(data = 0, nrow = num_row, ncol = 1)
  for (i in 1:num_features) {
  feature1 <- data_df[,i]
  mean_val <- mean(feature1)
  #print(mean_val)
  results_df[i,2] <- mean_val    
  }
  return(results_df)
}
 


#Read in quantile data
 
quantile_df <- read.csv("../data/all_quantile_df.csv")
#View(quantile_df)
 

 
ACPA_neg_df <- filter(quantile_df, quantile_df$Sample_Status == 'ACPA_negative')
#View(ACPA_neg_df)
ACPA_pos_df <- filter(quantile_df, quantile_df$Sample_Status == 'ACPA_positive')
#View(ACPA_pos_df)
control_df <- filter(quantile_df, quantile_df$Sample_Status == 'control')
#View(control_df)
#Number of samples ACPA-: 30, #ACPA+: 32, #Control: 30
 

#Gather the features names only
 
feature_names_only <- data.frame(colnames(ACPA_pos_df[,1:1622]))
#View(feature_names_only)
 

#Dataframe for ACPA_pos
 
ACPA_pos_means_df <- data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
colnames(ACPA_pos_means_df) <- cbind("Autoantibodies","ACPA_pos_Mean")
ACPA_pos_means_df[,1] <- feature_names_only
#View(ACPA_pos_means_df)
 

#Dataframe for ACPA_neg
 
ACPA_neg_means_df <- data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
colnames(ACPA_neg_means_df) <- cbind("Autoantibodies","ACPA_neg_Mean")
ACPA_neg_means_df[,1] <- feature_names_only
#View(ACPA_neg_means_df)
 

#Dataframe for control
 
control_means_df <- data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
colnames(control_means_df) <- cbind("Autoantibodies","control_Mean")
control_means_df[,1] <- feature_names_only
#View(control_means_df)
 

#Fill the ACPA_pos dataframe
 
feature1 <- matrix(data = 0, nrow = 32, ncol = 1)
for (i in 1:1622) {
  feature1 <- ACPA_pos_df[,i]
  mean_val <- mean(feature1)
  ACPA_pos_means_df[i,2] <- mean_val
}
#View(ACPA_pos_means_df)
 

#Find the feature means for ACPA_pos
 
num_features1 = ncol(ACPA_pos_df)-2
#num_features
ACPA_pos_means_df <- data.frame(matrix(data = 0, nrow = 1622, ncol = 2))
colnames(ACPA_pos_means_df) <- cbind("Autoantibodies","ACPA_pos_Mean")
ACPA_pos_means_df[,1] <- feature_names_only
ACPA_pos_means_df <- feature_means(ACPA_pos_means_df, ACPA_pos_df, num_features1)
#View(ACPA_pos_means_df)

 

#Now fill ACPA_neg and control means
 
ACPA_neg_means_df <- feature_means(ACPA_neg_means_df, ACPA_neg_df, num_features1)
#View(ACPA_neg_means_df)
control_means_df <- feature_means(control_means_df, control_df, num_features1)
#View(control_means_df)
 

#I need to combine the means to plot.
 
ACPA_pos_and_control_means <- cbind(ACPA_pos_means_df,control_means_df[,2])
colnames(ACPA_pos_and_control_means) <- cbind("Autoantibodies", "ACPA_pos_means", "control_means")
#View(ACPA_pos_and_control_means)
 

#Now do the same for ACPA_neg
#I need to combine the means of ACPA_pos and control to plot.
 
ACPA_neg_and_control_means <- cbind(ACPA_neg_means_df,control_means_df[,2])
colnames(ACPA_neg_and_control_means) <- cbind("Autoantibodies", "ACPA_neg_means", "control_means")
#View(ACPA_neg_and_control_means)
 

#Be sure to include all autoantibodies
 
all_RA_log_changes <- data.frame(ACPA_pos_and_control_means)
all_RA_log_changes$ACPA_neg_means <- ACPA_neg_and_control_means$ACPA_neg_means
all_RA_log_changes$ACPA_pos_log_change <- 0
all_RA_log_changes$ACPA_neg_log_change <- 0
#View(all_RA_log_changes)
 

#Fill in the log change values
 
for (i in 1:num_features1) {
  sero_mean <- all_RA_log_changes[i,2]
  control_mean <- all_RA_log_changes[i,3]
  seropos_log <- log2(sero_mean/control_mean)
  all_RA_log_changes[i,5] <- seropos_log
  
  seroneg_mean <- all_RA_log_changes[i,4]
  seroneg_log <- log2(seroneg_mean/control_mean)
  all_RA_log_changes[i,6] <- seroneg_log
} 
 

#add a column for the log change status, if the absolute value for both #is 
#greater than log2(1.5) then I will color it.
 
all_RA_log_changes$log_change_status <- "No significant change"
 

#Find significant changes higher in RA
 
for (i in 1:num_features1) {
  ACPA_pos_log <- all_RA_log_changes[i,5]
  ACPA_neg_log <- all_RA_log_changes[i,6]
  absolute_ACPA_pos <- abs(ACPA_pos_log)
  absolute_ACPA_neg <- abs(ACPA_neg_log)
  if(absolute_ACPA_pos > 0.5849625 | absolute_ACPA_neg > 0.5849625){
    all_RA_log_changes[i,7] <- "Higher in RA"
  }
}
 

#Plot
 
scatterplot1 <- ggplot(all_RA_log_changes, aes(x= ACPA_pos_log_change,y=ACPA_neg_log_change, color = log_change_status )) + geom_point()+ geom_abline(intercept = 0, slope = 1) +xlim(-2.5, 2.5)+ ylim(-2.5,2.5)+ xlab("log2(ACPA+/Control)")+ ylab("log2(ACPA-/Control") + coord_fixed() + scale_colour_manual(values = c("red4","black"))
scatterplot1
pdf("../output/scatter_plot1.pdf")
scatterplot1
dev.off()

#Ternary Plot
 
library("ggtern")
library("plotly")
library("readr")
library("dplyr")
library("tidyr")
 

#First I need to group all the means for the three groups
 
all_mean_values_df <- cbind(ACPA_pos_means_df,ACPA_neg_means_df[,2],control_means_df[,2])
colnames(all_mean_values_df) <- cbind("Autoantibody","ACPA_pos", "ACPA_neg", "Control")
#View(all_mean_values_df)
 

#Plot
 
ternary_plot1 <- ggtern(data = all_mean_values_df, aes(x=ACPA_neg, y = ACPA_pos, z=Control))+geom_point()+theme_bw()

pdf("../output/ternary_plot1.pdf")
ternary_plot1
dev.off()

#Plot with axis labels
 
ternary_plot2 <- ggtern(data = all_mean_values_df, aes(x=ACPA_neg, y = ACPA_pos, z=Control))+geom_point()+theme_bw()+theme_rgbw()

pdf("../output/ternary_plot2.pdf")
ternary_plot2
dev.off() 




