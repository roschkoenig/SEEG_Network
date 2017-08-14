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
### Prepare files for DCM analysis (SPM: MEEG object)
```
seeg_filemaker
```
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
