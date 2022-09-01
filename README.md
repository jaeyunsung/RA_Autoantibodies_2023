RA Autoantibodies 2022
===================================================
# Introduction
This repository contains the source code to reproduce the main data presented in Cunningham *et al.*, "Patients with ACPA-positive and ACPA-negative Rheumatoid Arthritis Show Different Serological Autoantibody Repertoires and Autoantibody Associations with Disease Activity". 

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

# Data used in this study

The raw autoantibody abundance (RFU) data and the quantile normalized data are all available in the 'data' directory.

# Description of scripts

The 'scripts' directory contains several R scripts for running the analysis pipelines for this project. Each script calls local data from the 'data' directory, and all results used to generate the main manuscript's display items are automatically put into the
'output' directory.

The scripts found in this repository are the following:

>differentially_abundant.R

This script takes in the raw autoantibody abundance data and transforms the RFU values using quantile normalization.
Next, the script finds statistically significant autoantibodies using a Mann-Whitney *U* test
and the Cliff's delta effect size. Significant features have a *P*-value less than
0.05 and an absolute Cliff's delta value greater than 0.33 (*i.e.*, moderate effect size). There 
are four grid plots that are made from running this script (Figures 2 and 3). Each grid plot shows the autoantibodies' abundances 
in an RA subgroup and controls.

>Figure1A_heatmap.R

This script takes the quantile normalized data and performs z-score transformation.
This transformed dataset is used to make a heatmap of all 1,622 autoantibody abundances for all samples.

>Figure1B_PCA.R

This script takes the quantile normalized data and first runs an ANOVA test to find features that
have statistically significant variance across study groups. If the *P*-value is less than 0.05, then the feature 
is considered as statistically significant. There are 26 significant features that are plotted in the
principal component analysis (PCA) plot.


>Figure1C_D_ternary_scatter_plots.R

This script uses the quantile normalized autoantibody abundances to build a ternary plot and a fold-change plot.
The fold-change plot is between ACPA+ RA/Control and ACPA– RA/Control.


>Figure5_bubble_plot.R

This script takes the quantile normalized data and the corresponding RA patients' CDAI. The script finds
the Spearman correlation between the individual autoantibody abundances and CDAI scores. ACPA+ RA, ACPA–
RA, and all RA patients (*i.e.*, ACPA+ RA and ACPA– RA) are the three groups for which we calculate Spearman correlations.
The significant autoantibodies (*P*-value < 0.01 & |rho| > 0.4) are plotted in the bubble plot.

# Running the scripts

Example of how to run a script in the terminal command line:

> Rscript differentially_abundant.R
