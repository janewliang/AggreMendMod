library(abind)

# Helper functions
source("../diagnostics_functions.R")

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
short_cancers = cancer_map[PanelPRO:::MODELPARAMS$PanPRO22$CANCERS]

# Number of relatives with > 1 cancer in each family
mult_cancer_rels = sapply(results, function(x) {
  sum(rowSums(x$fam[paste0("isAff", short_cancers)]) > 1)
})

###############################################################################

# Pull out the proband information from each simulated family
# Mutation status
mut_df = data.frame(t(sapply(results, function(x){
  unlist(x$fam[x$fam$isProband==1,PanelPRO:::MODELPARAMS$PanPRO22$GENES])
})))
mut_df$ALLGENES = ifelse(rowSums(mut_df) > 0, 1, 0)

# Pull out the probabilities for each model
prob_df = lapply(1:nrow(results[[1]]$probs), function(i){
  data.frame(t(sapply(results, function(x){x$probs[i,]})))
})
names(prob_df) = c("PP22_All", "PP22_First", "Agg")

save(prob_df, mut_df, mult_cancer_rels, file="results/output/all_dfs.rData")

###############################################################################

# Get diagnostics
all_out = list(
  AllFamilies = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES, x$ALLGENES)
  }), 
  MultipleCancers1 = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[mult_cancer_rels > 0], 
                        x$ALLGENES[mult_cancer_rels > 0])
  }), 
  MultipleCancers0 = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[mult_cancer_rels == 0], 
                        x$ALLGENES[mult_cancer_rels == 0])
  })
)

diagnostics = lapply(all_out, function(x) {
  sapply(x, function(y){y[[1]]})
})
carrier_probs = lapply(all_out, function(x) {
  lapply(x, function(y){y[2:length(y)]})
})

# Save output
save(diagnostics, mult_cancer_rels, 
     file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
