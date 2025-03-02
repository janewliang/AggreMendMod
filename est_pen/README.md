Code to estimate aggregate penetrances. 

Files: 
- `estimate_pp21.R`: Code to estimate penetrances aggregating over 21 genes and 17 cancers, based on simulated data. 

Sub-directories: 
- `pp21/`: Simulate families and estimate aggregate penetrances under the assumptions of Fam3PRO-21, a 21-gene and 17-cancer model. 

Each sub-directory contains the following files: 
- `fam.R`: Code to simulate 1000 families at a time.  
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for simulating families, for usage on a high performance computing cluster. 
- `estimate.sh`: Shell script for estimating penetrances, for usage on a high performance computing cluster. 
