
 
library(dplyr)
library(reshape2)
library(ggpubr)
library(Hmisc)
library(preprocessCore)
library(effsize)
 

#Functions for analysis
 
analysis_part1_plot = function(temp_df, plot_title){
#plot all samples as box and whisker
  
    temp_plot = ggboxplot(temp_df, x = "Case_ID", y = "value",color = 'black',
                       title = plot_title, short.panel.labs = FALSE,
                       outlier.shape=NA, xlab = "Samples", ylab = "Auto-antibody levels", 
                       fill = "Sample_Status") + scale_y_continuous(limits=c(-10, 1000))+theme( axis.text.x = element_blank())
    return (temp_plot)
}

make_result_table = function(matrix_1, matrix_2, num_features, feature_list){
#This will do Wilcox and CliffDelta simultaneously
    #Define result table
    result_table = as.data.frame(matrix(data = 0, nrow = num_features, ncol = 3))
    result_table[,1] = feature_list
    colnames(result_table) = cbind("Autoantibody", "Nominal_pvalue", "Cliff_delta")
    #Do Wilcoxon
    for(i in 1:num_features){
        return_val = wilcox.test(matrix_1[,i], matrix_2[,i])[3]
        result_table[i,2] = return_val
    }
    #Do CliffDelta
    for(i in 1:num_features){
        cliff_val = cliff.delta(d = matrix_1[,i], f = matrix_2[,i])[1]
        result_table[i,3] = cliff_val
    }
    #Sort by Pval
    result_table = result_table[order(result_table$Nominal_pvalue),]
    return (result_table)
}

get_sig_feature_list = function(wilcox_cliff_result_table){
#Find feature list that has p < 0.05 & abs(cliffdelta) > 0.33
    temp_list = c()
    num_features = nrow(wilcox_cliff_result_table)
    for(i in 1:num_features){
        feature = wilcox_cliff_result_table$Autoantibody[i]
        pval = wilcox_cliff_result_table$Nominal_pvalue[i]
        cliff_delta = wilcox_cliff_result_table$Cliff_delta[i]
        cliff_delta = abs(cliff_delta)
        if (pval < 0.05){
            if (cliff_delta > 0.33){
                temp_list = c(temp_list, feature)
            }
        }
    }
    return (temp_list)  
}

autoantibody_boxplots <- function(sig_features_df){
  #take a dataframe of features and make box plots comparing both groups
  temp_df <- melt(sig_features_df)
  box_plot <- ggplot(data = temp_df, aes(x = variable, y =value, color = Sample_Status)) + geom_boxplot() + facet_wrap(~variable, scales="free", ncol = 5) +stat_compare_means(label = "p.format", method=  "wilcox.test", label.y.npc =0.9) +theme(strip.background = element_blank(),strip.text.x = element_blank()) + labs( x = "Auto-antibody")+ scale_colour_manual(values = c("brown1", "steelblue"))#+ geom_jitter()
  return(box_plot)
}

get_sigfeatures = function(num_sigfigs, sigfeature_list, data_df, results_df){
  #Gather the significant features into a new dataframe
  for (i in 1:num_sigfigs) {
    feature1 <- sigfeature_list[i,1]
    feature_to_add <- data_df[feature1]
    results_df[,i] <- feature_to_add
  }
  return(results_df)
}

 

#Read raw data in


normal_data_df <- read.csv("../data/normal_data_rows_as_samples.csv")
#View(normal_data_df)
 

#Feature names only
 
num_cols <- ncol(normal_data_df)
feature_names_df <- colnames(normal_data_df[2:num_cols])
#View(feature_name_df)
 

#Class labels for samples
 
class_labels_df = read.csv("../data/class_labels_with_sample_ID_no_spaces.csv")
colnames(class_labels_df) <- cbind("Case_ID", "Sample_Status")
#View(class_labels_df)
 

#Get the number of columns for later
 
#num_cols
total_samples <- nrow(normal_data_df)
#total_samples
 

#Combine Case IDs and sample status
 
normal_data_rows_as_samples_to_plot = merge(normal_data_df, class_labels_df, by='Case_ID', all=TRUE)
#View(normal_data_rows_as_samples_to_plot)
dim(normal_data_rows_as_samples_to_plot)
 

#Find sample medians for all data
 
normal_data_rows_as_samples_to_plot$Column_Medians <- apply(normal_data_rows_as_samples_to_plot[,2:num_cols], 1, median)
#View(normal_data_rows_as_samples_to_plot)
 

#Sort by median
 
normal_data_rows_as_samples_to_plot <- data.frame(normal_data_rows_as_samples_to_plot[order(normal_data_rows_as_samples_to_plot$Column_Medians, decreasing = FALSE),])
#View(normal_data_rows_as_samples_to_plot)


#Plot all samples by median
 
melt_df1 <- melt(normal_data_rows_as_samples_to_plot) #Note: First time calling reshape library
all_samples_plot <- analysis_part1_plot(melt_df1,"Sample distribution before Normalization")
all_samples_plot
 

#Remove top 5
 
reduced_normal_data_df <- normal_data_rows_as_samples_to_plot[1:(nrow(normal_data_rows_as_samples_to_plot)-5),1:(num_cols+1)]
reduced_df <- melt(reduced_normal_data_df)
reduced_samples_plot <- analysis_part1_plot(reduced_df,"Sample distribution after reduction")
reduced_samples_plot
 

#Two patients were on rituximab treatment so we will remove them
 
only_92_samples_df <- subset(reduced_normal_data_df, Case_ID != "Case08" & Case_ID != "Case59") 
#View(only_92_samples_df)


#Do quantile normalization
 
remove_list = c("Sample_Status", "Case_ID")
before_quantile_df = only_92_samples_df[,which(colnames(only_92_samples_df) %nin% remove_list)] 
#This will replace LINE 172 - LINE 200

transpose_before_quantile_df = t(before_quantile_df)
after_quantile_df = normalize.quantiles(transpose_before_quantile_df)
quantile_92_df = data.frame(t(after_quantile_df))

colnames(quantile_92_df) <- feature_names_df
#View(quantile_92_df)

quantile_92_df$Case_ID = only_92_samples_df$Case_ID
quantile_92_df$Sample_Status = only_92_samples_df$Sample_Status

quantiled_data_df = melt(quantile_92_df)
quantile_plot <- analysis_part1_plot(quantiled_data_df, "Sample distribution after quantile normalization")
quantile_plot
 

#Export the data for later use
 
#write.csv(quantile_df, "../RA_Autoantibodies_2022/QN_sengenics_data_with_class.csv", row.names=FALSE, quote=FALSE)
#write.csv(quantile_92_df, "/Users/Kevin/Desktop/all_quantile_df.csv", row.names=FALSE, quote=FALSE)

#Seperate the three groups
 

ACPA_neg_df <- filter(quantile_92_df, quantile_92_df$Sample_Status == 'ACPA_negative')
#View(ACPA_neg_df)
num_ACPA_neg <- nrow(ACPA_neg_df)

ACPA_pos_df <- filter(quantile_92_df, quantile_92_df$Sample_Status == 'ACPA_positive')
#View(ACPA_pos_df)
num_ACPA_pos <- nrow(ACPA_pos_df)

control_df <- filter(quantile_92_df, quantile_92_df$Sample_Status == 'control')
#View(control_df)
num_control <- nrow(control_df)
 

#Some variables to use
 
num_features <- ncol(quantile_92_df) - 2
#num_features
feature_list <- feature_names_df
#feature_list
 

#ACPA_pos or ACPA_neg is the first group. 
#ALl positive Cliff's Delta values will indicate higher abundance in RA. 
#All negative Cliff's Delta values will indicate higher abundance in control. 
 
posVcont_result_table <- make_result_table(ACPA_pos_df, control_df, num_features, feature_list)
#View(posVcont_result_table)
negVcont_result_table <- make_result_table(ACPA_neg_df, control_df, num_features, feature_list)
#View(negVcont_result_table)
#save the output
#write.csv(posVcont_result_table, "../RA_Autoantibodies_2022/posVcont_result_table.csv", row.names=FALSE, quote=FALSE)
#write.csv(negVcont_result_table, "../RA_Autoantibodies_2022/negVcont_result_table", row.names=FALSE, quote=FALSE)
 

#Select features significant between ACPA_pos and control
 
higher_in_ACPApos <- subset(posVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta > 0.33)
#View(higher_in_ACPApos)

higher_in_control_thanACPApos <- subset(posVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta < -0.33)
#View(higher_in_control_thanACPApos)
 

#Gather the features higher in ACPA_pos than in control
 
#how many significant features for ACPA_pos
num_sigfigs_pos <- nrow(higher_in_ACPApos)
#num_sigfigs_pos
ACPA_pos_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_pos, ncol = num_sigfigs_pos))
colnames(ACPA_pos_sigfeats_df1) <- higher_in_ACPApos[,1]
ACPA_pos_sigfeats_df1 <- get_sigfeatures(num_sigfigs_pos, higher_in_ACPApos, ACPA_pos_df, ACPA_pos_sigfeats_df1)
ACPA_pos_sigfeats_df1$Sample_Status <- ACPA_pos_df$Sample_Status
#View(ACPA_pos_sigfeats_df1)
 

#Gather the features in control that are lower than ACPA_pos
 
num_sigfigs_controlPos <- nrow(higher_in_ACPApos)
ACPA_pos_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_controlPos))
colnames(ACPA_pos_sigfeats_df2) <- higher_in_ACPApos[,1]
ACPA_pos_sigfeats_df2 <- get_sigfeatures(num_sigfigs_pos, higher_in_ACPApos, control_df, ACPA_pos_sigfeats_df2)
ACPA_pos_sigfeats_df2$Sample_Status <- control_df$Sample_Status
#View(ACPA_pos_sigfeats_df2)
 

#Combine and plot
 
high_in_ACPA_pos_than_control_df <- rbind(ACPA_pos_sigfeats_df1 ,ACPA_pos_sigfeats_df2 )
#View(high_in_ACPA_pos_than_control_df)
ACPA_posControl_plot <- autoantibody_boxplots(high_in_ACPA_pos_than_control_df)

pdf("../output/higher_in_ACPA_pos_than_control.pdf")
ACPA_posControl_plot
dev.off() 

#Now get the features that are higher in control compared to 
#ACPA_pos from the ACPA_pos data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlpos <- nrow(higher_in_control_thanACPApos)
control_higher_than_pos_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_pos, ncol = num_sigfigs_Controlpos))
colnames(control_higher_than_pos_sigfeats_df1) <- higher_in_control_thanACPApos[,1]
control_higher_than_pos_sigfeats_df1 <- get_sigfeatures(num_sigfigs_Controlpos, higher_in_control_thanACPApos, ACPA_pos_df, control_higher_than_pos_sigfeats_df1)
control_higher_than_pos_sigfeats_df1$Sample_Status <- ACPA_pos_df$Sample_Status
#View(control_higher_than_pos_sigfeats_df1)
 

#Now get the features that are higher in control compared to
#ACPA_pos from the control data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlpos <- nrow(higher_in_control_thanACPApos)
control_higher_than_pos_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_Controlpos))
colnames(control_higher_than_pos_sigfeats_df2) <- higher_in_control_thanACPApos[,1]
control_higher_than_pos_sigfeats_df2 <- get_sigfeatures(num_sigfigs_Controlpos, higher_in_control_thanACPApos, control_df, control_higher_than_pos_sigfeats_df2)
control_higher_than_pos_sigfeats_df2$Sample_Status <- control_df$Sample_Status
#View(control_higher_than_pos_sigfeats_df2)
 

#Combine and plot
 
higher_in_control_than_in_ACPA_pos_df <- rbind(control_higher_than_pos_sigfeats_df1, control_higher_than_pos_sigfeats_df2 )
#View(higher_in_control_than_in_ACPA_pos_df)
Control_ACPA_pos_plot <- autoantibody_boxplots(higher_in_control_than_in_ACPA_pos_df)

pdf("../output/higher_in_control_than_ACPA_pos.pdf")
Control_ACPA_pos_plot
dev.off()  

#Select features significant between ACPA_neg and control
 
higher_in_ACPAneg <- subset(negVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta > 0.33)
#View(higher_in_ACPAneg)

higher_in_control_thanACPAneg <- subset(negVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta < -0.33)
#View(higher_in_control_thanACPAneg)
 

#Gather the features higher in ACPA_neg than in control
 
#how many significant features for ACPA_pos
num_sigfigs_neg <- nrow(higher_in_ACPAneg)
ACPA_neg_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_neg, ncol = num_sigfigs_neg))
colnames(ACPA_neg_sigfeats_df1) <- higher_in_ACPAneg[,1]
ACPA_neg_sigfeats_df1 <- get_sigfeatures(num_sigfigs_neg, higher_in_ACPAneg, ACPA_neg_df, ACPA_neg_sigfeats_df1)
ACPA_neg_sigfeats_df1$Sample_Status <- ACPA_neg_df$Sample_Status
#View(ACPA_neg_sigfeats_df1)
 

#Gather the same features in control
 
num_sigfigs_neg <- nrow(higher_in_ACPAneg)
ACPA_neg_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_neg))
colnames(ACPA_neg_sigfeats_df2) <- higher_in_ACPAneg[,1]
ACPA_neg_sigfeats_df2 <- get_sigfeatures(num_sigfigs_neg, higher_in_ACPAneg, control_df, ACPA_neg_sigfeats_df2)
ACPA_neg_sigfeats_df2$Sample_Status <- control_df$Sample_Status
#View(ACPA_neg_sigfeats_df2)
 

#Combine and plot
 
higher_in_ACPA_neg_than_in_control_df <- rbind(ACPA_neg_sigfeats_df1, ACPA_neg_sigfeats_df2)
#View(higher_in_ACPA_neg_than_in_control_df)
ACPA_neg_control_plot <- autoantibody_boxplots(higher_in_ACPA_neg_than_in_control_df)

pdf("../output/higher_in_ACPA_neg_than_control.pdf")
ACPA_neg_control_plot
dev.off()

#Now get the features that are higher in control compared to 
#ACPA_neg from the ACPA_neg data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlneg <- nrow(higher_in_control_thanACPAneg)
control_higher_than_neg_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_neg, ncol = num_sigfigs_Controlneg))
colnames(control_higher_than_neg_sigfeats_df1) <- higher_in_control_thanACPAneg[,1]
control_higher_than_neg_sigfeats_df1 <- get_sigfeatures(num_sigfigs_Controlneg, higher_in_control_thanACPAneg, ACPA_neg_df, control_higher_than_neg_sigfeats_df1)
control_higher_than_neg_sigfeats_df1$Sample_Status <- ACPA_neg_df$Sample_Status
#View(control_higher_than_neg_sigfeats_df1)
 

#Now get the features that are higher in control compared to
#ACPA_pos from the control data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlneg <- nrow(higher_in_control_thanACPAneg)
control_higher_than_neg_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_Controlneg))
colnames(control_higher_than_neg_sigfeats_df2) <- higher_in_control_thanACPAneg[,1]
control_higher_than_neg_sigfeats_df2 <- get_sigfeatures(num_sigfigs_Controlneg, higher_in_control_thanACPAneg, control_df, control_higher_than_neg_sigfeats_df2)
control_higher_than_neg_sigfeats_df2$Sample_Status <- control_df$Sample_Status
#View(control_higher_than_neg_sigfeats_df2)
 

#Combine and plot
 
higher_in_control_than_in_ACPA_neg_df <- rbind(control_higher_than_neg_sigfeats_df1, control_higher_than_neg_sigfeats_df2 )
#View(higher_in_control_than_in_ACPA_neg_df)
Control_ACPA_neg_plot <- autoantibody_boxplots(higher_in_control_than_in_ACPA_neg_df)
Control_ACPA_neg_plot
 
pdf("../output/higher_in_control_than_ACPA_neg.pdf")
Control_ACPA_neg_plot
dev.off()


