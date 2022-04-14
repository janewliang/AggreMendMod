# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("../run_model_functions.R")

# Load HCP families (1362)
load("hcp_families.rData")

###############################################################################

# Split the families for parallelization
idx = split(1:length(hcp_families_PP), 
            ceiling(seq_along(1:length(hcp_families_PP)) / 
                      ceiling(length(hcp_families_PP)/681)))


# Run models for comparison on families
fam_output = lapply(idx[[s]], function(i){
  # Extract family
  fam_PP = hcp_families_PP[[i]]
  
  # Add extra columns
  # Get short cancer names
  cancer_map = PanelPRO:::CANCER_NAME_MAP$short
  names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long
  short_cancers = c(cancer_map[PanelPRO:::MODELPARAMS$PanPRO22$CANCERS], 
                    "OtherCancers")
  fam_PP$isAffAllCancers = ifelse(rowSums(fam_PP[,paste0("isAff", short_cancers)], 
                                          na.rm = TRUE) > 0, 1, 0)
  fam_PP$AgeAllCancers = apply(fam_PP[,paste0("Age", short_cancers)], 1, 
                               function(x) {
                                 temp = min(x, na.rm = TRUE)
                                 temp[is.infinite(temp)] = NA
                                 return(temp)
                               })
  
  # Run models
  out = try(run_models(fam_PP), TRUE)
  
  # Return carrier probabilities
  return(out)
})


# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
