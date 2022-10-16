# DiauxicShift

Untargeted metabolomics for functional discovery of regulatory genes in the context of the diauxic shift.
This repository contains code and intermediate data to be able to recreate the work if needed. Supplementary data can be found in /data and /results.

Last update: 2022-10-16

## Project information

Semi-automated batch cultivations of the yeast S. cerevisiae (BY4741, MATa his3&#916;1 leu2&#916;0 met15&#916;0 ura3&#916;0) and 10 single deletion mutants (&#916;ygro67c, &#916;tda1, &#916;rme1, &#916;rts3, &#916;pcl1, &#916;oca1, &#916;gal11, &#916;dld3, &#916;mek1, &#916;faa1) sampled in the fermentative phase and respiratory phase respectively. Untargeted LC-MS metabolomics was done on the processed samples and subsequently analyzed in order to gain increased understanding of the diauxic shift and the functionality of its regulators. 

Supplementary material, such as processed data and strainwise statistics can be found inside this repository in the data-folder.
Raw data in the form of mzML-files can be found at Zenodo: 10.5281/zenodo.7105589

## FBA Simulations

Simulations done to select the strains of interest was made using the framework supplied in Coutant et al. (https://doi.org/10.1073/pnas.1900548116). Output found in data/simulations.

## Data processing

Analysis makes use of the output provided by MSDial (v4.7) as a basis for processing and analysis. The scripts should be ran in the following order:

1. scripts/processing_step0.py (Reads MSDial outputs and performs basic cleaning and curation)
2. scripts/processing_step1.R (Reads cleaned output and performs signal drift correction and batch correction)
3. scripts/processing_step2.R (Basic curation and data formatting)
4. scripts/processing_step3.R (Merges MS-runs)
5. scripts/processing_MAnalyst.R (Prepares data for downstream MetaboAnalyst analysis, outputs ft.csv (data file) and st.csv (metadata file) for the phase and strain analysis)

## Analysis

Basic analysis made using the MetaboAnalyst online-client (https://www.metaboanalyst.ca/)

Analysis using FELLA and visualizations performed in the following scripts:

1. scripts/figure1A.R (Heatmap of metabolite profiles)
2. scripts/figure1BC_data.R (FELLA-based analysis of metaboanalyst output)
3. scripts/figure1BC_plot.R (Visualization of enrichment)
4. scripts/figure2B.R (Volcano plot of differentially expressed metabolites)
5. scripts/figure2C.R (FELLA-based analysis of the diauxic shift)
6. scripts/figure4.R (Visualization of spearman correlation of simulations and statistically significant differences)
