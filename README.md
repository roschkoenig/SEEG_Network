# Network modelling of stereotactic EEG data 
*work in progress* | [About this research stream](http://www.dynamic-brains.com/networks-in-epilepsy-surgery/) | [ICTALS2017 Poster](https://doi.org/10.6084/m9.figshare.5311198.v1)

This repository contains the code required to reproduce a network based approach to SEEG data analysis that aims to simulate epilepsy surgery *in silico* to help improve surgical planning in the future. This project is **work in progress** and has not yet been validated for prospective clinical use.
The strategy demonstrated here rests on three broad approahces
1. Using nonnegative matrix factorisation (an idea stolen from [Chai et al. 2017](http://doi.org/10.1162/NETN_a_00001)) to easily visualise **network edge dynamics** during seizures in a patient undergoing presurgical evaluation for epilepsy surgery using stereotactic EEG. 
2. Using hierarchical dynamic causal modelling (parametric empirical Bayesian DCM, which rests on [Friston et al. 2016](https://doi.org/10.1016/j.neuroimage.2015.11.015)) to identify **biophysical network changes** that best explain seizure onset across a number of different seizure types within the same patient. 
3. Using simulations of fully fitted neural mass models (the output of the DCM analysis) to quantify the effects of **in silico node resection** on whole network activity

The code runs on [Matlab](https://uk.mathworks.com/) (tested with 2016b) and requires the following freely available packages to run
* [Statistical Parametric Mapping, SPM12](http://www.fil.ion.ucl.ac.uk/spm/) - The neuroimaging analysis suite that contains the DCM tools

## Table of Contents
1. [Custom routines](#custom-routines-included-in-this-repository)
    1. [Visualise Network](#visualise-dynamic-network-patterns)
    1. [Nonnegative Matrix Decomposition](#run-nonnegative-matrix-decomposition-to-identify-subgraphs)
    1. [Subnetwork Channels](#identifying-a-subnetwork-of-representative-channels)
    1. [MEEG Object Filemaker](#prepare-files-for-dcm-analysis-spm-meeg-object)
    1. [Grand Mean DCM](#run-grand-mean-dcm-to-specify-prior-parameter-estimates)
    1. [DCMs of seizures](#run-dcm-on-regions-of-interest-identified-from-nonnegative-matrix-decomposition)
    1. [Parametric Empirical Bayes](#run-hierarchical-parametric-empirical-model-to-identify-seizure-related-model-changes)
    1. [Simulated Surgery](#simulate-surgical-interventions-and-quantify-network-effects)
2. [SPM Modifications]
3. [Other Helper Functions]


## Custom routines included in this repository

The repository includes a number of custom routines that when run sequentially should reproduce our analysis. 

### Visualise dynamic network patterns
```
seeg_networkvis
```
In order to first visualise the network dynamics during seizures in this dataset, this routine plots some of the activity. The first section plots individual channels for the duraction of the seizure activity. 
<img src="https://user-images.githubusercontent.com/12950773/29360266-15c9d3f8-827a-11e7-9823-7fdb68354fe8.png">

The routine then employs a sliding window approach to estimate the changing output spectral distribution at each channel over time (note that the sliding window in its current implementation will take a [**long** time](https://github.com/roschkoenig/SEEG_Network/issues/2)). 

<img src="https://user-images.githubusercontent.com/12950773/29378229-32a5129e-82b6-11e7-90e9-3d9b14250b56.png">

The routine also plots example adjacency matrices derived from windowed coherence and correlation estimates for a pre-seizure and seizure segment in the data. 

<img src="https://user-images.githubusercontent.com/12950773/29378230-32b40b3c-82b6-11e7-8ebf-86a20fa79741.png" width="500">

To track the full temporal resolution of these changes, we can vectorise the adjacency matrix (using the matlab inbuilt function `tril` to preserve just the unique edges) and plot the resultant vectors for each time window. Note that for better visibility, this routine orders the edges according to their overall mean correlation/coherence. 

<img src="https://user-images.githubusercontent.com/12950773/29378595-994477f0-82b7-11e7-9590-08b5b5abcb08.png">





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
1. Multichannel segmented into sliding window, so that each 'trial' in the MEEG-object corresponds to a separate *30s* window in the original, concatenated data. This will be used to perform [window-of-interest DCM analysis](#run-dcm-on-regions-of-interest-identified-from-nonnegative-matrix-decomposition). 
2. Single channel concatenated data, so that each 'trial' in the MEEG-object corresponds to a separate channel from the 7-channel selection chosen for DCM analysis. This will be used in the [Grand Mean DCM](#run-grand-mean-dcm-to-specify-prior-parameter-estimates) analysis. 

The screenshot below shows the first (multichannel, windowed) MEEG object when seen in SPM (open in matlab by calling `spm eeg`, then select M/EEG from the 'Display' drop down menu. 

<img src="https://user-images.githubusercontent.com/12950773/29273228-6c3ab7b4-80fb-11e7-9050-5de90f9daed8.PNG" alt="Expression Peaks" width="400">

### Run grand mean DCM to specify prior parameter estimates
```
seeg_gmdcm
```
Prior parameters for the biophysical models used in SPM are based on cortical columns in healthy human cortex, and are specifically designed to cope well with ERP experiments and particular steady state applications. Here we need to adapt the priors to allow the inversion of abnormal dynamics as seen in epilepsy. 

In order to find suitable sections of parameter space, we therefore employ a staged approach to the model inversion
1. Manually set a small selection of parameters to values that are compatible with the spectral output seen in the data 
2. Use this parametersiation to invert continuous data from individual channels, yielding a single, 'steady-state' neural mass model parameterisation that reproduces the type of spectral output seen in this particular channel
3. Use the channel-specific parameterisations as priors for further DCM analysis of multichannel networks identifying changes 

The first two points are addressed in this routine. Running`spm_induced_optimise` allows visualisation of individual parameters' effects on te spectral output of individual model nodes. Based on this, a single parameter (the time constant *T(2)*) was altered. With the thus adapted priors, this routine runs an inversion for each channel individually. The empirical estimates for the individual channels is then saved for use as priors in the next steps of the analysis. 

### Run DCM on regions of interest identified from [Nonnegative Matrix Decomposition](#identifying-a-subnetwork-of-representative-channels)
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
