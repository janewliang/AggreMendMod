library(PanelPRO)
library(abind)

# Helper functions
source("../diagnostics_functions.R")

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

# 443 sets of family results
S = 443

# Putting all of the results together
all_probs = NULL
for (s in 1:S){
  load(paste0("results/output/fam", s, ".rData"))
  all_probs = c(all_probs, fam_output)
}

# 7 failed families out of 443
failed_families = sapply(all_probs, function(x){length(x)==0 || class(x)=="try-error"})
failed_errs = all_probs[failed_families]
names(failed_errs) = which(failed_families)
all_probs[failed_families] = NULL

# Identify the panel with 21 genes tested
panel_id_21 = which(sapply(panel_list, length) == 21)

# Only keep families with 21 genes tested (591)
mgh_families_PP = mgh_families_PP[proband_muts$PanelID == panel_id_21]
mgh_families_PP[failed_families] = NULL

###############################################################################

# Pull out the proband mutation information from each family
proband_muts = proband_muts[proband_muts$PanelID == panel_id_21,]
mut_df = proband_muts[!failed_families,genes]
mut_df[,"ALLGENES"] = ifelse(rowSums(mut_df, na.rm = TRUE) >= 1, 1, 0)

# Pull out the probabilities for each model
prob_df = lapply(1:nrow(all_probs[[1]]), function(i) {
  data.frame(t(sapply(all_probs, function(x){x[i,]})))
})
names(prob_df) = c("PP21", "Agg")

save(prob_df, mut_df, file="results/output/all_dfs.rData")

###############################################################################

# Get diagnostics
all_out = lapply(prob_df, function(x) {
  get_all_diagnostics(mut_df$ALLGENES, x$ALLGENES)
})

diagnostics = sapply(all_out, function(x) {x[[1]]})
carrier_probs = lapply(all_out, function(x) {x[2:length(x)]})

# Save output
save(diagnostics, failed_errs, failed_families, 
     file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
