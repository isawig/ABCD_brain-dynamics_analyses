# ABCD_brain-dynamics_analyses
Codes for all the analyses implemented in the "Functional connectome harmonics and dynamic connectivity maps of the preadolescent brain" manuscript.

*A_prepare-for-connectomics/* folder contains script to run before connectome analysis. 
-  Input data are preprocessed resting-state mean timeseries, with ROIs as columns and timepoints as rows.
-  Other inputs required are: .txt/.csv file with subjects to keep for the analysis (in this case with framewise displacement, FD < 0.4).
-  The output are subjects mean timeseries, with ROIs as rows and timepoints as columns; with first 100 cortical ROIs (Schaefer100) and 14 subcortical means (to match ENIGMA Toolbox requirements).

*B_LEiDA-analysis/* folder contains main steps to run LEiDA from the whole sample, derive useful metrics, and validate network results with reference resting-state networks (RSNs).

-  *block-A_run-LEiDA-on-HPC/* contains scripts to run LEiDA on CSC high-performance computer (HPC) environment. Input data are subjects mean timeseries; and k-means clustering results, as well as fractional occupancies and dwell times estimates, are given as output.
-  *block-B_derive-metrics/* contains scripts to derive Markov Chain transition probabilities starting from fractional occupancies (FO) and dwell times (DT) estimated in block A.
-  *block-C_network-validation/* 
