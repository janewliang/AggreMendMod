Code to run simulation studies.  

Files: 
- `diagnostics_functions.R`: Functions for obtaining diagnostic metrics.
- `diagnostics_boot.R`: Code to obtain diagnostic metrics from 10 bootstrap samples.
- `combine_diagnostics_boot.R`: Code to load bootstrapped diagnostic metrics returned from separate cluster jobs.

Sub-directories: 
- `pp21/`: Simulate and evaluate families under the assumptions of Fam3PRO-21, a 21-gene and 17-cancer model. 
- `cancer_misreporting/`: Simulate and evaluate families under the assumptions of Fam3PRO-21 but with cancer type misreporting. 

Each sub-directory contains the following files: 
- `functions.R`: Functions that run models and extract carrier probability results. 
- `fam.R`: Code to simulate 1000 families at a time and run models. 
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for running models on the simulated families, for usage on a high performance computing cluster. 
- `diagnostics.R`: Code to obtain diagnostic metrics from evaluating models on simulated families.
- `diagnostics.sh`: Shell script for obtaining diagnostic metrics, for usage on a high performance computing cluster. 
- `rscript_diagnostics_boot_job`and `submitJobs_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrap samples, for usage on a high performance computing cluster. 
- `analysis.R`: Code to generate figures and tables. 
