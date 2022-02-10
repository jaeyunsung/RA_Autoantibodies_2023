RA Autoantibodies 2022

This repository includes several R scripts as well R markdown scripts to run
for the analysis for this project. The input data in included in the repository
and no editing is required to run the scripts. Each script has local data
that is read in and all output for the figures is automatically put into the
output directory.

Each R script has a correponding R markdown file that does the same thing.
The scripts included are the following:

differentially_abundant_final.R
This script take in the raw data and quantile normalizes the samples first.
Next, the script finds significant autoantibodies using a Mann-Whitney U test
and Cliff's delta effect size. Significant features have a P-value less than
0.05 and an absolute Cliff's delta value greater than 0.33. There are four plots
that are made from running this script. Each shows the RA abundance and the control
abundance for comparison.

Figure1_heatmap_final.R
This script takes the quantile data and does a z-score transformation on the data.
Then it makes a heatmap of all 1,622 features for all samples.

Figure1_PCA_final.R
This script takes the quantile data and first runs an ANOVA test to find features that
are statistically significant. If the P-value is less than 0.05 then the feature is
considered significant. There are 26 significant features that are then plotted as a
PCA plot.

Figure1_scatter_final.R
The Clinical Disease Activity score is a commonly used RA measurement of disease
severity. In this script, the CDAI scores are used with the quantile normalized data
to find Spearman Correlations. The autoantibody with the highest Spearman Correlation is
for the CISH gene.

Installation
There are a few R libraries that need to be installed to run these scripts.
They are:
library(dplyr)
library(reshape2)
library(ggpubr)
library(Hmisc)
library(preprocessCore)
library(effsize)
library(pheatmap)
library(FactoMineR)
library(stats)
library(factoextra)
library(ggfortify)

NOTE: You must make a directory called "output" in the same directory as all the scripts and files.
This is where the figures will be saved.

Here is an example of running the scripts:

Rscript differentially_abundant_final.r
