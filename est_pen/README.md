Code to estimate aggregate penetrances. 

Files: 
- `estimate_pp22.R`: Code to estimate penetrances aggregating over 22 genes and 17 cancers, based on simulated data. 
- `estimate_oc.R`: Code to estimate penetrances aggregating over 22 genes and 17 cancers, plus "any other cancer", based on simulated data. 

Sub-directories: 
- `pp22/`: Simulate families and estimate aggregate penetrances under the assumptions of PanelPRO-22, a 22-gene and 17-cancer model. 
- `oc*/`, where `*` is replaced by `1.2`, `1.5`, `1.8`, `2.1`, `2.4`, `2.7`, or `3`: Simulate families and estimate aggregate penetrances under the assumptions of PanelPRO-22, a 22-gene and 17-cancer model, plus "any other cancer". `*` corresponds to the multiplying factor used to compute the carrier and noncarrier penetrances for "any other cancer". 

Each sub-directory contains the following files: 
- `fam.R`: Code to simulate 1000 families at a time.  
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for simulating families, for usage on a high performance computing cluster. 
- `estimate.sh`: Shell script for estimating penetrances, for usage on a high performance computing cluster. 

In addition, the `oc*/` sub-directories contain the following file: 
- `get_database.R`: Code to set up a database with "any other cancer" that is compatible with the `PanelPRO` R package. 