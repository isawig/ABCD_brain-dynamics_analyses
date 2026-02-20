# ABCD_brain-dynamics_analyses
Codes for all the analyses implemented in the "Functional connectome harmonics and dynamic connectivity maps of the preadolescent brain" manuscript.

## A_prepare-for-connectomics
*A_prepare-for-connectomics/* folder contains script to run before connectome analysis. 
-  Input data are preprocessed resting-state mean timeseries, with ROIs as columns and timepoints as rows.
-  Other inputs required are: .txt/.csv file with subjects to keep for the analysis (in this case with framewise displacement, FD < 0.4).
-  The output are subjects mean timeseries, with ROIs as rows and timepoints as columns; with first 100 cortical ROIs (Schaefer100) and 14 subcortical means (to match ENIGMA Toolbox requirements).

## B_LEiDA-analysis
*B_LEiDA-analysis/* folder contains main steps to run LEiDA from the whole sample, derive useful metrics, and validate network results with reference resting-state networks (RSNs).

-  *block-A_run-LEiDA-on-HPC/* contains scripts to run LEiDA on CSC high-performance computer (HPC) environment. Input data are subjects mean timeseries; and k-means clustering results, as well as fractional occupancies and dwell times estimates, are given as output.
-  *block-B_derive-metrics/* contains scripts to derive Markov Chain transition probabilities starting from fractional occupancies (FO) and dwell times (DT) estimated in block A.
-  *block-C_network-validation/* estimates correlations between LEiDA-derived centroids and RSNs binary vectors, comparing only regions whose phase aligns to the main eigenvector (positive centroids value) and assigns each centroid to the best matching network. Then results are visualized on the cortical surface with ENIGMA Toolbox, color-coded according to the assigned RSN, using for correlations a significance threshold of r > 0.4. It takes as input LEiDA centroids and FOs and returns correlations and networks assignment, as well as pyramid visualization.

## C_FCH-analysis
*C_FCH-analysis/* folder contains main step to run FCH analysis from the whole sample, derive power and energy metrics, and compare first low-order harmonics (low-frequency associated ones) with seven adults reference gradients.

-  *block-A_run_FCH_pipeline/* allows to run FCH pipeline from subjects mean timeseries, and to derive functional harmonics as well as cortical projections, useful to estimate following metrics.
-  *block-B_derive-metrics/* contains code to estimate power and energy from harmonics and eigenvectors.
-  *block-C_PLSR_with_LEiDA/* takes as input LEiDA centroids for given K, and FCH ordered and randomized harmonics, and fits them in a PLSR model. It returns as output the MSE, percentage of explained variance, and other model estimates.
-  *block-D_NeuroMap_validation/* compares low-order FCHs to seven reference adults gradients via correlation analysis, and returns as output the estimated correlations between each pair of brain networks. It takes as input one harmonic at a time and compares it to all gradients. Then results are visualized with ENIGMA. 

## D_statistical-analysis
*D_statistical-analysis/* is where harmonization procedure is performed, and connectome-derived metrics are used for statistical analyses.









