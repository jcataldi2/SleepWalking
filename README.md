# SleepWalking

The script "Perform_glmm_CE_EEG.R" performs the glmm (generalized linear mixed model) to predict CE (conscious experience) vs. NE (no experience) or vs. CEWR (conscious experience without recall) for each EEG channel/voxel.
based on the power spectral density (PSD) calculated on EEG signal in different timeframes and 2 frequencies (Delta and Beta), using the identity of the subject as random factor.
INPUT: .xlsx file where each line is a channel/voxel. The columns contain the info of the subject identity, trial identity, CE/NE/CEWR info (CE = 2, CEWR = 1, NE = 0), a columns containing the PSD for each frequency and each time frame.

OUTPUT: a R. data containing the statistical values (beta values, z-values and p values) per channel/voxel per frequency.

To plot the statistical results (z-values and p-values):
- channel level: use the "topoplot" eeglab function
- source level (voxels): use the "plotInflatedMap_SL_new" function.
