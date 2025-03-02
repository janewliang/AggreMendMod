# Extracting the task and job IDs
s = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Helper functions
source("../run_model_functions.R")

# Load MGH families
load("../../../Clinical_validation/mgh/data/pp21_other_cancers/mgh_families.rData")

# Genes and cancers
genes =  c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", 
           "CDH1", "CDK4", "CDKN2A", "CHEK2", "EPCAM", 
           "MLH1", "MSH2", "MSH6", "NBN", "PALB2", "PMS2", 
           "PTEN", "RAD51C", "RAD51D", "STK11", "TP53")
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", "Gastric", 
            "Kidney", "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", 
            "Pancreas", "Prostate", "Small Intestine", "Soft Tissue Sarcoma", 
            "Thyroid", "Urinary Bladder", "Hepatobiliary")

###############################################################################

# Only keep families with 21 genes tested (443)
panel_id_21 = which(sapply(panel_list, length) == 21)
mgh_families_PP = mgh_families_PP[proband_muts$PanelID == panel_id_21]

# Split the families for parallelization
idx = split(1:length(mgh_families_PP), 
            ceiling(seq_along(1:length(mgh_families_PP)) / 
                      ceiling(length(mgh_families_PP)/443)))


# Run models for comparison on families
fam_output = lapply(idx[[s]], function(i){
  # Extract family
  fam_PP = mgh_families_PP[[i]]
  
  # Add extra columns
  # Get short cancer names
  cancer_map = PanelPRO:::CANCER_NAME_MAP$short
  names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long
  short_cancers = c(cancer_map[cancers], "OtherCancers")
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
