Two Matlab tools used by various projects to compute cluster-based permutation statistics on time series data.


eegck_clusterstats.m opeartes on D-dimensional data and compares a specific test statistics between the actual data
and randomized data. The data is assumed to be continous in D-dimensions.
Different criteria for cluster-definition and first- and second-level thresholds can be provided. 
Requires fieldtrip or SPM to be installed.


eegck_clusterstats_eeg.m opeartes on EEG data. An additional input providing the fieldtrip electrode positions
is required, and clusters are computed between neighbouring electrode-time points. 
