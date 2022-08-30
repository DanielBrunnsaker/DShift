# DiauxicShift

Untargeted metabolomics for functional discovery of regulatory genes in the context of the diauxic shift.
This repository contains code and intermediate data to be able to recreate the work.

Last update: 2022-04-14

## Project information

Semi-automated batch cultivations of the yeast S. cerevisiae (BY4741, MATa his3&#916;1 leu2&#916;0 met15&#916;0 ura3&#916;0) and 10 single deletion mutants (&#916;ygro67c, &#916;tda1, &#916;rme1, &#916;rts3, &#916;pcl1, &#916;oca1, &#916;gal11, &#916;dld3, &#916;mek1, &#916;faa1) sampled in the fermentative phase and respiratory phase respectively. Untargeted LC-MS metabolomics was done on the processed samples and subsequently analyzed in order to gain increased understanding of the diauxic shift and the functionality of its regulators.

## FBA Simulations

Simulations done to select the strains of interest was made using the framework supplied in Coutant et al. (https://doi.org/10.1073/pnas.1900548116). Output found in data/simulations.

## Data processing

Analysis makes use of the output provided by MSDial (v4.7) as a basis for processing and analysis. The scripts should be ran in the following order:

1. processing_step1.R
2. processing_step2.R
3. processing_step3.R
4. processing_MAnalyst.R
