library(PanelPRO)
library(abind)

# Helper functions
source("../diagnostics_functions.R")

# Load MGH families
load("mgh_families.rData")

###############################################################################

# 591 sets of family results
S = 591

# Putting all of the results together
all_probs = NULL
for (s in 1:S){
  load(paste0("results/output/fam", s, ".rData"))
  all_probs = c(all_probs, fam_output)
}

# 7 failed families out of 591
failed_families = sapply(all_probs, function(x){length(x)==0 || class(x)=="try-error"})
failed_errs = all_probs[failed_families]
names(failed_errs) = which(failed_families)
all_probs[failed_families] = NULL

# Identify the panel with 22 genes tested
panel_id_22 = which(sapply(panel_list, length) == 22)

###############################################################################

# Only keep families with 22 genes tested (591)
mgh_families_PP = mgh_families_PP[proband_muts$PanelID == panel_id_22]
mgh_families_PP[failed_families] = NULL

# Get short cancer names
cancer_map = PanelPRO:::CANCER_NAME_MAP$short
names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long
short_cancers = c(cancer_map[PanelPRO:::MODELPARAMS$PanPRO22$CANCERS], 
                  "OtherCancers")

# Number of relatives with > 1 cancer in each family
mult_cancer_rels = sapply(mgh_families_PP, function(fam) {
  sum(rowSums(fam[paste0("isAff", short_cancers)]) > 1)
})

# Proband ancestry
ancestry = sapply(mgh_families_PP, function(fam) {
  fam$Ancestry[fam$isProband==1]
})

###############################################################################

# Pull out the proband mutation information from each family
proband_muts = proband_muts[proband_muts$PanelID == panel_id_22,]
mut_df = proband_muts[!failed_families,PanelPRO:::MODELPARAMS$PanPRO22$GENES]
mut_df[,"ALLGENES"] = ifelse(rowSums(mut_df, na.rm = TRUE) >= 1, 1, 0)

# Pull out the probabilities for each model
prob_df = lapply(1:nrow(all_probs[[1]]), function(i) {
  data.frame(t(sapply(all_probs, function(x){x[i,]})))
})
# Multiplying factors for OtherCancers carrier penetrance
multiplying_factors = seq(1.2, 3, by = 0.3)
names(prob_df) = c("PP22_RM", "PP22_All", "PP22_First", "Agg_no_oc", 
                   paste0("Agg_", multiplying_factors))

save(prob_df, mut_df, mult_cancer_rels, ancestry, 
     file="results/output/all_dfs.rData")

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
  }), 
  nonAJ = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[ancestry == "nonAJ"], 
                        x$ALLGENES[ancestry == "nonAJ"])
  }), 
  AJ_Italian = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[ancestry != "nonAJ"], 
                        x$ALLGENES[ancestry != "nonAJ"])
  }), 
  AJ = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[ancestry == "AJ"], 
                        x$ALLGENES[ancestry == "AJ"])
  }), 
  Italian = lapply(prob_df, function(x) {
    get_all_diagnostics(mut_df$ALLGENES[ancestry == "Italian"], 
                        x$ALLGENES[ancestry == "Italian"])
  })
)

diagnostics = lapply(all_out, function(x) {
  sapply(x, function(y){y[[1]]})
})
carrier_probs = lapply(all_out, function(x) {
  lapply(x, function(y){y[2:length(y)]})
})

# Save output
save(diagnostics, failed_errs, failed_families, mult_cancer_rels, 
     file="results/diagnostics/diagnostics.rData")
save(carrier_probs, file="results/diagnostics/carrier_probs.rData")
