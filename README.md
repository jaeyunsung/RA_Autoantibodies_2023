RA Autoantibodies 2022
===================================================
# Introduction
This repository contains the source code to reproduce the main data presented in "Patients with ACPA-positive and ACPA-negative Rheumatoid Arthritis Show Different Serological Autoantibody Repertoires and Autoantibody Associations with Disease Activity" by Cunningham *et al.* (manuscript submitted for publication). 

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

# Description of data used in this study

The 'data' directory contains the following input files for the analysis pipelines of this project:

>all_quantile_df.csv

Quantile normalized autoantibody abundance data.

>patient_info_table.csv

Clinical Disease Activity Index (CDAI) scores of RA patients. Of note, only CDAI scores above 0 (n = 56) are shown.

>quantile_for_CDAI_scatterplots.csv

Quantile normalized autoantibody abundance data for RA patients whose CDAI scores are above 0 (n = 56).

>raw_92_samples.csv

Raw autoantibody abundance (RFU) data.

# Description of scripts

The 'scripts' directory contains several R scripts for running the analysis pipelines of this project. Each script calls local data from the 'data' directory. All results used to generate the main text's display items are automatically exported into the 'output' directory.

The scripts found in this repository are the following:

>differentially_abundant.R

This script takes in the raw autoantibody abundance data and transforms the RFU values using quantile normalization.
Next, differentially abundant autoantibodies are found using the Mann-Whitney *U* test
and the Cliff's delta effect size. Statistically significant features have a *P*-value less than
0.05 and a Cliff's delta magnitude greater than 0.33 (*i.e.*, moderate effect size). There 
are four grid plots that are made from running this script (Figures 2 and 3). Each grid plot shows the autoantibodies' abundances 
in an RA subgroup and controls.


>Figure1A_heatmap.R

This script takes the quantile normalized data and performs z-score transformation.
The transformed dataset matrix is then used to make a heatmap displaying 1,622 autoantibody abundances for the 92 profiles.


>Figure1B_PCA.R

This script takes the quantile normalized data and first runs an ANOVA test to find features (i.e., autoantibodies) that
have statistically significant variance across study groups. If the *P*-value is less than 0.05, then the feature 
is considered to be statistically significant. Only the significant features are used to project autoantibody profiles onto the
principal component analysis (PCA) ordination plot.


>Figure1C_D_ternary_scatter_plots.R

This script uses the quantile normalized autoantibody abundances to build a ternary plot and a fold-change scatter-plot.
The scatter-plot compares fold-changes of ACPA+ RA vs. Controls and of ACPA– RA vs. Controls.


>Figure5_bubble_plot.R

This script takes the quantile normalized data and the corresponding RA patients' CDAI scores. The Spearman correlation (rho) between the individual autoantibody abundances and CDAI scores are found.  Of note, only CDAI scores above 0 (n = 56) are used to find correlations with autoantibody abundances. ACPA+ RA, ACPA–
RA, and all RA patients (*i.e.*, combining both ACPA+ RA and ACPA– RA into a single group) are the three groups for which we calculate Spearman correlations.
Only the significant autoantibodies (*P*-value < 0.01 & |*rho*| > 0.4) are shown in the bubble plot.

# Running the scripts

Example of how to run a script in the terminal command line:

> Rscript differentially_abundant.R
