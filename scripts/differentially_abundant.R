library(dplyr)
library(reshape2)
library(ggpubr)
library(Hmisc)
library(preprocessCore)
library(effsize)

#Functions for analysis
 
analysis_part1_plot = function(temp_df, plot_title){
#plot all samples as box and whisker
  
    temp_plot = ggboxplot(temp_df, x = "Sample_ID", y = "value",color = 'black',
                       title = plot_title, short.panel.labs = FALSE,
                       outlier.shape=NA, xlab = "Samples", ylab = "Auto-antibody levels", 
                       fill = "Study_group") + scale_y_continuous(limits=c(-10, 1000))+theme( axis.text.x = element_blank())
    return (temp_plot)
}

make_result_table = function(matrix_1, matrix_2, num_features, feature_list){
#“This will perform a Mann-Whitney U test and identify the Cliff’s delta effect size”
    #Define result table
    result_table = as.data.frame(matrix(data = 0, nrow = num_features, ncol = 3))
    result_table[,1] = feature_list
    colnames(result_table) = cbind("Autoantibody", "Nominal_pvalue", "Cliff_delta")
    #Do Mann-Whitney U test
    for(i in 1:num_features){
        return_val = wilcox.test(matrix_1[,i], matrix_2[,i])[3]
        result_table[i,2] = return_val
    }
    #Do Cliff's delta
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
  box_plot <- ggplot(data = temp_df, aes(x = variable, y =value, color = Study_group)) + geom_boxplot() + facet_wrap(~variable, scales="free", ncol = 5) +stat_compare_means(label = "p.format", method=  "wilcox.test", label.y.npc =0.9) +theme(strip.background = element_blank(),strip.text.x = element_blank()) + labs( x = "Auto-antibody")+ scale_colour_manual(values = c("brown1", "steelblue"))#+ geom_jitter()
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

#Read in data
print("read in the raw data.")
raw_df <- read.csv("../data/raw_92_samples.csv")


#Feature names only
 
num_cols <- ncol(raw_df)-2
#print(num_cols)
feature_names_df <- colnames(raw_df[1:num_cols])
#print(feature_names_df)


#Quantile normalization
remove_list = c("Study_group", "Sample_ID")
before_quantile_df = raw_df[,which(colnames(raw_df) %nin% remove_list)] 
#print(before_quantile_df)

transpose_before_quantile_df = t(before_quantile_df)
after_quantile_df = normalize.quantiles(transpose_before_quantile_df)
quantile_df = data.frame(t(after_quantile_df))
#print(quantile_df)

colnames(quantile_df) <- feature_names_df
quantile_df$Sample_ID = raw_df$Sample_ID
quantile_df$Study_group = raw_df$Study_group
write.csv(quantile_df, "../data/all_quantile_df.csv", row.names = FALSE)
quantiled_data_df = melt(quantile_df)
quantile_plot <- analysis_part1_plot(quantiled_data_df, "Sample distribution after quantile normalization")
print(quantile_plot)

#Seperate the three groups
 

ACPA_neg_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_negative_RA')
#View(ACPA_neg_df)
num_ACPA_neg <- nrow(ACPA_neg_df)
print("number of ACPA negative: ")
print(num_ACPA_neg)

ACPA_pos_df <- filter(quantile_df, quantile_df$Study_group == 'ACPA_positive_RA')
#View(ACPA_pos_df)
num_ACPA_pos <- nrow(ACPA_pos_df)
print("number of ACPA positive: ")
print(num_ACPA_pos)

control_df <- filter(quantile_df, quantile_df$Study_group == 'control')
#View(control_df)
num_control <- nrow(control_df)
print("number of control: ")
print(num_control)

num_features <- ncol(quantile_df) - 2
#print(num_features)
feature_list <- feature_names_df
#print(feature_list)

#ACPA_pos or ACPA_neg is the first group. 
#ALl positive Cliff's Delta values will indicate higher abundance in RA. 
#All negative Cliff's Delta values will indicate higher abundance in control. 
posVcont_result_table <- make_result_table(ACPA_pos_df, control_df, num_features, feature_list)
#print(posVcont_result_table)

negVcont_result_table <- make_result_table(ACPA_neg_df, control_df, num_features, feature_list)

higher_in_ACPApos <- subset(posVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta > 0.33)
#print(higher_in_ACPApos)
higher_in_control_thanACPApos <- subset(posVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta < -0.33)
#print(higher_in_control_thanACPApos)

#Gather the features higher in ACPA_pos than in control
 
#how many significant features for ACPA_pos
num_sigfigs_pos <- nrow(higher_in_ACPApos)
#num_sigfigs_pos
ACPA_pos_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_pos, ncol = num_sigfigs_pos))
colnames(ACPA_pos_sigfeats_df1) <- higher_in_ACPApos[,1]
ACPA_pos_sigfeats_df1 <- get_sigfeatures(num_sigfigs_pos, higher_in_ACPApos, ACPA_pos_df, ACPA_pos_sigfeats_df1)
ACPA_pos_sigfeats_df1$Study_group <- ACPA_pos_df$Study_group
#print(ACPA_pos_sigfeats_df1)


#Gather the features in control that are lower than ACPA_pos
 
num_sigfigs_controlPos <- nrow(higher_in_ACPApos)
ACPA_pos_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_controlPos))
colnames(ACPA_pos_sigfeats_df2) <- higher_in_ACPApos[,1]
ACPA_pos_sigfeats_df2 <- get_sigfeatures(num_sigfigs_pos, higher_in_ACPApos, control_df, ACPA_pos_sigfeats_df2)
ACPA_pos_sigfeats_df2$Study_group <- control_df$Study_group
#print(ACPA_pos_sigfeats_df2)


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
control_higher_than_pos_sigfeats_df1$Study_group <- ACPA_pos_df$Study_group
#print(control_higher_than_pos_sigfeats_df1)


#Now get the features that are higher in control compared to
#ACPA_pos from the control data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlpos <- nrow(higher_in_control_thanACPApos)
control_higher_than_pos_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_Controlpos))
colnames(control_higher_than_pos_sigfeats_df2) <- higher_in_control_thanACPApos[,1]
control_higher_than_pos_sigfeats_df2 <- get_sigfeatures(num_sigfigs_Controlpos, higher_in_control_thanACPApos, control_df, control_higher_than_pos_sigfeats_df2)
control_higher_than_pos_sigfeats_df2$Study_group <- control_df$Study_group
#View(control_higher_than_pos_sigfeats_df2)



#Combine and plot
 
higher_in_control_than_in_ACPA_pos_df <- rbind(control_higher_than_pos_sigfeats_df1, control_higher_than_pos_sigfeats_df2 )
#View(higher_in_control_than_in_ACPA_pos_df)
Control_ACPA_pos_plot <- autoantibody_boxplots(higher_in_control_than_in_ACPA_pos_df)

pdf("../output/lower_in_ACPA_pos_than_control.pdf")
Control_ACPA_pos_plot
dev.off()  



#Select features significant between ACPA_neg and control
 
higher_in_ACPAneg <- subset(negVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta > 0.33)
#print(higher_in_ACPAneg)

higher_in_control_thanACPAneg <- subset(negVcont_result_table, Nominal_pvalue < 0.05 & Cliff_delta < -0.33)
#print(higher_in_control_thanACPAneg)


#Gather the features higher in ACPA_neg than in control
 
#how many significant features for ACPA_pos
num_sigfigs_neg <- nrow(higher_in_ACPAneg)
ACPA_neg_sigfeats_df1 <- data.frame(matrix(data = 0, nrow = num_ACPA_neg, ncol = num_sigfigs_neg))
colnames(ACPA_neg_sigfeats_df1) <- higher_in_ACPAneg[,1]
ACPA_neg_sigfeats_df1 <- get_sigfeatures(num_sigfigs_neg, higher_in_ACPAneg, ACPA_neg_df, ACPA_neg_sigfeats_df1)
ACPA_neg_sigfeats_df1$Study_group <- ACPA_neg_df$Study_group
#View(ACPA_neg_sigfeats_df1)


#Gather the same features in control
 
num_sigfigs_neg <- nrow(higher_in_ACPAneg)
ACPA_neg_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_neg))
colnames(ACPA_neg_sigfeats_df2) <- higher_in_ACPAneg[,1]
ACPA_neg_sigfeats_df2 <- get_sigfeatures(num_sigfigs_neg, higher_in_ACPAneg, control_df, ACPA_neg_sigfeats_df2)
ACPA_neg_sigfeats_df2$Study_group <- control_df$Study_group
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
control_higher_than_neg_sigfeats_df1$Study_group <- ACPA_neg_df$Study_group
#View(control_higher_than_neg_sigfeats_df1)


#Now get the features that are higher in control compared to
#ACPA_pos from the control data
 
#how many significant features higher in control than in ACPA_pos
num_sigfigs_Controlneg <- nrow(higher_in_control_thanACPAneg)
control_higher_than_neg_sigfeats_df2 <- data.frame(matrix(data = 0, nrow = num_control, ncol = num_sigfigs_Controlneg))
colnames(control_higher_than_neg_sigfeats_df2) <- higher_in_control_thanACPAneg[,1]
control_higher_than_neg_sigfeats_df2 <- get_sigfeatures(num_sigfigs_Controlneg, higher_in_control_thanACPAneg, control_df, control_higher_than_neg_sigfeats_df2)
control_higher_than_neg_sigfeats_df2$Study_group <- control_df$Study_group
#View(control_higher_than_neg_sigfeats_df2)

#Combine and plot
 
higher_in_control_than_in_ACPA_neg_df <- rbind(control_higher_than_neg_sigfeats_df1, control_higher_than_neg_sigfeats_df2 )
#View(higher_in_control_than_in_ACPA_neg_df)
Control_ACPA_neg_plot <- autoantibody_boxplots(higher_in_control_than_in_ACPA_neg_df)
Control_ACPA_neg_plot
 
pdf("../output/lower_in_ACPA_neg_than_control.pdf")
Control_ACPA_neg_plot
dev.off()
file.remove("./Rplots.pdf")
