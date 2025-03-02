library(abind)

# Helper functions
source("../diagnostics_functions.R")

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

# 1000 sets of family results
K = 1000

# Putting all of the family simulation results together
results = NULL
for (k in 1:K){
  load(paste0("results/output/fam", k, ".rData"))
  results = c(results, fam_output)
}

###############################################################################

# Get short cancer names
cancer_map = PanelPRO:::CANCER_NAME_MAP$short
names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long
short_cancers = cancer_map[cancers]

###############################################################################

# Pull out the proband information from each simulated family
# Mutation status
mut_df = data.frame(t(sapply(results, function(x){
  unlist(x$fam[x$fam$isProband==1,genes])
})))
mut_df$ALLGENES = ifelse(rowSums(mut_df) > 0, 1, 0)

# Pull out the probabilities for each model
prob_df = unlist(lapply(1:length(c("10", "30", "50")), function(j) {
  lapply(1:nrow(results[[1]]$probs[[1]]), function(i){
    data.frame(t(sapply(results, function(x){x$probs[[j]][i,]})))
  })
}), recursive = FALSE)
names(prob_df) = paste0(rep(c("PP21", "Agg"), 2), "_", 
                        rep(c(10, 30, 50), each = 2))

save(prob_df, mut_df, file="results/output/all_dfs.rData")

###############################################################################

# Get diagnostics
all_out = lapply(prob_df, function(x) {
  get_all_diagnostics(mut_df$ALLGENES, x$ALLGENES)
})

diagnostics = sapply(all_out, function(x) {x[[1]]})
carrier_probs = lapply(all_out, function(x) {x[2:length(x)]})

# Save output
save(diagnostics, file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
