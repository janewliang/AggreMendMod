# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("functions.R")

# Load simulation functions
files.sources = paste0("../../simulate_families/", 
                       list.files("../../simulate_families/"))
sapply(files.sources, source)

# Load CGN family structure
load("../../../CGN/CGN.fam.structure.RData")

###############################################################################

# 1000 family simulations
nsim = 1000

# Run models for comparison on simulated families
fam_output = lapply(1:nsim, function(i){
  
  # Set random seed
  set.seed(s*nsim+i)
  # Sample from CGN families
  rn = sample(1:nrow(cgn.fam.structure), 1)
  
  # Randomly sample number of males and females in each branch of the family.
  # Paternal aunt, paternal uncles
  nSibsPatern = c(cgn.fam.structure[rn,3], cgn.fam.structure[rn,2]) 
  # Maternal aunts, maternal uncles
  nSibsMatern = c(cgn.fam.structure[rn,5], cgn.fam.structure[rn,4]) 
  # Sisters and brothers
  nSibs = c(cgn.fam.structure[rn,7], cgn.fam.structure[rn,6]) 
  # We make the assumption that the number of sons and daughters for the 
  # proband and all siblings, is the same. Nieces and nephews of the proband 
  # are not sampled separately.
  nGrandchild = c(cgn.fam.structure[rn,9], cgn.fam.structure[rn,8]) 
  
  # Simulate family
  fam_PP = sim.runSimFam(nSibsPatern, nSibsMatern, nSibs, nGrandchild, 
                         new_db, PanelPRO:::MODELPARAMS$PanPRO22$GENES, 
                         PanelPRO:::MODELPARAMS$PanPRO22$CANCERS, 
                         includeGeno = TRUE, affTime = TRUE, ageMin = 1)
  
  # Add extra columns
  fam_PP$isAffAllCancers = fam_PP$isAffAny
  fam_PP$AgeAllCancers = fam_PP$AgeAny
  
  # Swap cancer types for 10% of affected relatives
  fam_PP_10 = swapCancerTypes(fam_PP, 0.1)
  # Swap cancer types for 30% of affected relatives
  fam_PP_30 = swapCancerTypes(fam_PP, 0.3)
  # Swap cancer types for 50% of affected relatives
  fam_PP_50 = swapCancerTypes(fam_PP, 0.5)
  
  # Run models
  out = lapply(1:3, function(i) {
    run_models(list(fam_PP_10, fam_PP_30, fam_PP_50)[[i]])
  })
  names(out) = c(10, 30, 50)
  
  return(list(fam = fam_PP, probs = out))
})

# Save results
save(fam_output, file = paste0("results/output/fam", s, ".rData"))
