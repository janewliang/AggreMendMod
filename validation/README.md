Code to validate models on data from the USC-Stanford Hereditary Cancer Panel (HCP) Testing Study<sup>[1](#myfootnote1)</sup> and the Center for Cancer Risk Assessment at Massachusetts General Hospital (MGH). 

Files: 
- `run_model_functions.R`: Functions that run models and extract carrier probability results. 
- `diagnostics_functions.R`: Functions for obtaining diagnostic metrics.
- `combine_diagnostics_boot.R`: Code to load bootstrapped diagnostic metrics returned from separate cluster jobs.

Sub-directories: 
- `hcp/`: Validate the aggregate model on the HCP cohort. Assumes that pre-processed pedigrees are saved as a list in `hcp_families.rData`. 
- `mgh/`: Validate the aggregate model on the MGH cohort. Assumes that pre-processed pedigrees are saved as a list in `mgh_families.rData`. 

Each sub-directory contains the following files: 
- `fam.R`: Code to run models on the validation cohort. 
- `rscript_fam.job` and`submitJobs_fam.sh`: Shell scripts for running models on validation data, for usage on a high performance computing cluster. 
- `diagnostics.R`: Code to obtain diagnostic metrics from evaluating models on validation data.
- `diagnostics.sh`: Shell script for obtaining diagnostic metrics, for usage on a high performance computing cluster. 
- `diagnostics_boot.R`: Code to obtain diagnostic metrics from 10 bootstrap samples.
- `rscript_diagnostics_boot_job`and `submitJobs_diagnostics_boot.sh`: Shell scripts for obtaining diagnostic metrics from 1000 bootstrap samples, for usage on a high performance computing cluster. 
- `analysis.R`: Code to generate figures and tables. 


---

<a name="myfootnote1">1</a>. Idos, G. E., Kurian, A. W., Ricker, C., Sturgeon, D., Culver, J. O., Kingham, K. E., ... & Levonian, P. (2019). Multicenter prospective cohort study of the diagnostic yield and patient experience of multiplex gene panel testing for hereditary cancer risk. JCO Precision Oncology, 3, 1-12.
