# Network modelling of stereotactic EEG data

This repository contains the code required to reproduce a network based approach to SEEG data analysis that aims to simulate epilepsy surgery *in silico* to help improve surgical planning in the future. This project is **work in progress** and has not yet been validated for prospective clinical use.
The strategy demonstrated here rests on three broad approahces
1. Using nonnegative matrix factorisation (an idea stolen from [Chai et al. 2017](http://doi.org/10.1162/NETN_a_00001)) to easily visualise **network edge dynamics** during seizures in a patient undergoing presurgical evaluation for epilepsy surgery using stereotactic EEG. 
2. Using hierarchical dynamic causal modelling (parametric empirical Bayesian DCM, which rests on [Friston et al. 2016](https://doi.org/10.1016/j.neuroimage.2015.11.015)) to identify **biophysical network changes** that best explain seizure onset across a number of different seizure types within the same patient. 
3. Using simulations of fully fitted neural mass models (the output of the DCM analysis) to quantify the effects of **in silico node resection** on whole network activity

The code runs on [Matlab](https://uk.mathworks.com/) (tested with 2016b) and requires the following freely available packages to run
* [Statistical Parametric Mapping, SPM12](http://www.fil.ion.ucl.ac.uk/spm/) - The neuroimaging analysis suite that contains the DCM tools

## Custom routines included in this repository

The repository includes a number of custom routines that when run sequentially should reproduce our analysis. 

### Visualise dynamic network patterns
```
seeg_networkvis
```
### Run nonnegative matrix decomposition to identify subgraphs
```
seeg_decomp
```
### Identifying a subnetwork of representative channels
```
seeg_components
```
In order to keep the DCM analysis at a computationally tractable level of complexity, we need to select a subsection of channels that can be modelled. There are several ways to do this, and this solution here presents a simple approach that aims to (1) select original data channels so that the time series remain interpretable as raw EEG data, and (2) make the selection data-driven and automatic.

For this we use the same matrix-decomposition algorithm as above, only this time we use it on a symmetrical weighted adjacency matrix derived from coherence across all data points. This decomposition will give us two (very similar) *n by k* matrices, where *n* is the number of channels, and *k* is the number of components. Here we set the components to be identical to the number of SEEG electrodes in the original recording (i.e. 5), which produces the following results:

<img src="https://user-images.githubusercontent.com/12950773/29272047-480f66aa-80f6-11e7-8c77-15139af9481b.PNG" alt="Components" width="500">

We can plot these components on top of each other and identify regional maxima in the expression of specific components. These correspond to individual channels, and based on these peaks we select a small number of channels that are included in the subsequent DCM analysis. Note that this approach yields 7 peaks from the 5 components, which are included for further analysis in the DCM.  

<img src="https://user-images.githubusercontent.com/12950773/29272048-48269f46-80f6-11e7-8439-bafb1b749d74.PNG" alt="Expression Peaks" width="800">

### Prepare files for DCM analysis (SPM: MEEG object)
```
seeg_filemaker
```
SPM-based analyses use electrophysiological data in a specific matlab object category called the MEEG object. This routine creates two MEEG objects for use in the subsequent DCM analysis, which contains the data in two different ways
1. Multichannel segmented into sliding window, so that each 'trial' in the MEEG-object corresponds to a separate *30s* window in the original, concatenated data
2. Single channel concatenated data, so that each 'trial' in the MEEG-object corresponds to a separate channel from the 7-channel selection chosen for DCM analysis. This will be used in the [Grand Mean DCM](#Run-grand-mean-DCM-to-specify-prior-parameter-estimates) analysis. 

### Run grand mean DCM to specify prior parameter estimates
```
seeg_gmdcm
```
### Run DCM on regions of interest identified from [Nonnegative Matrix Decomposition](#Identifying-a-subnetwork-of-representative-channels)
```
seeg_dcm
```
### Run hierarchical parametric empirical model to identify seizure-related model changes
```
seeg_peb
```
### Simulate surgical interventions and quantify network effects
```
seeg_surgery
```
