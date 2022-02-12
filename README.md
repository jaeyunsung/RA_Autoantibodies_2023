RA Autoantibodies 2022
===================================================
# Introduction
This repository contains the source code to reproduce the main data presented in Cunningham et al., "Patients with ACPA-positive and ACPA-negative Rheumatoid Arthritis Show Differences in Circulating Autoantibody Repertoires". Manuscript in preparation.

# Description

The 'scripts' directory contains several R scripts for running the analysis pipelines for this project. Each script calls local data from the 'data' directory, and all results used to generate the main manuscript's display items are automatically put into the
'output' directory.

The scripts found in this repository are the following:

>differentially_abundant_final.R

This script take in the raw data and first transforms the data via quantile normalization.
Next, the script finds significant autoantibodies using a Mann-Whitney U test
and Cliff's delta effect size. Significant features have a P-value less than
0.05 and an absolute Cliff's delta value greater than 0.33. There are four plots
that are made from running this script. Each shows the RA abundance and the control
abundance for comparison.

>Figure1_heatmap_final.R

This script takes the quantile normalized data and performs z-score transformation.
This transformed dataset is used to make a heatmap of all 1,622 features for all samples.

>Figure1_PCA_final.R

This script takes the quantile normalized data and first runs an ANOVA test to find features that
are statistically significant. If the P-value is less than 0.05 then the feature is
considered significant. There are 26 significant features that are then plotted in the
principal components analysis (PCA) plot.

>Figure1_scatter_final.R

The Clinical Disease Activity (CDAI) score is a commonly used RA measurement of disease
activity. In this script, the CDAI scores are used with the quantile normalized data
to find Spearman correlation coefficients. The autoantibody with the highest Spearman correlation is
for the CISH gene.

# Data

The raw data and the quantile normalized data are all available in the 'data' directory.
There is also a text file named "removed_samples_info.txt" explaining why certain autoantibody samples (corresponding to study participants) were removed.

# Installation

The following R libraries must be pre-installed to run the scripts:

```
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
```

# Running the scripts

Here is an example in command line of how to run a script in the Linux Shell (terminal):

> Rscript differentially_abundant_final.r
