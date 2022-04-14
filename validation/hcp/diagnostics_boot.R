library(abind)

# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("../diagnostics_functions.R")

# Load results
load("results/output/all_dfs.rData")

###############################################################################

# Bootstrap diagnostics
R = 10 # Number of bootstrap replicates
n = nrow(mut_df) # Number of probands/families

# Get diagnostics from bootstrapped samples
diagnostics_boot = lapply(1:R, function(i){
  # Set random seed
  set.seed(s*R+i)
  
  # Sample observations with replacement
  idx = sample(n, n, replace=TRUE)
  
  # Bootstrap data
  boot_mut_df = mut_df[idx,]
  boot_prob_df = lapply(prob_df, function(x){x[idx,]})
  boot_mult_cancer_rels = mult_cancer_rels[idx]
  boot_race = race[idx]
  
  # Get bootstrap diagnostics
  all_out = list(
    AllFamilies = sapply(boot_prob_df, function(x) {
      get_all_diagnostics(boot_mut_df$ALLGENES, x$ALLGENES, 
                          return_probs=FALSE)
    }), 
    MultipleCancers1 = sapply(boot_prob_df, function(x) {
      get_all_diagnostics(boot_mut_df$ALLGENES[boot_mult_cancer_rels > 0], 
                          x$ALLGENES[boot_mult_cancer_rels > 0], 
                          return_probs=FALSE)
    }), 
    MultipleCancers0 = sapply(boot_prob_df, function(x) {
      get_all_diagnostics(boot_mut_df$ALLGENES[boot_mult_cancer_rels == 0], 
                          x$ALLGENES[boot_mult_cancer_rels == 0], 
                          return_probs=FALSE)
    }), 
    WNH = sapply(boot_prob_df, function(x) {
      get_all_diagnostics(boot_mut_df$ALLGENES[boot_race == "WNH"], 
                          x$ALLGENES[boot_race == "WNH"], 
                          return_probs=FALSE)
    }), 
    Hispanic = sapply(boot_prob_df, function(x) {
      get_all_diagnostics(boot_mut_df$ALLGENES[boot_race == "Hispanic"], 
                          x$ALLGENES[boot_race == "Hispanic"], 
                          return_probs=FALSE)
    })
  )
  return(all_out)
})

save(diagnostics_boot, file = paste0("results/diagnostics/diagnostics_boot", 
                                     s, ".rData"))
