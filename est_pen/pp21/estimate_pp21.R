library(PanelPRO)
library(abind)

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

# Pull out the proband information from each simulated family
# Mutation status
mut_df = data.frame(t(sapply(results, function(fam){
  unlist(fam[fam$isProband==1,genes])
})))
mut_df$ALLGENES = ifelse(rowSums(mut_df) > 0, 1, 0) 

# Cancer affection time without censoring
afftime_df = data.frame(t(sapply(results, function(fam){
  c(unlist(fam[fam$isProband==1, grep("Time*", names(fam))]), 
    unlist(fam[fam$isProband==1, c("Sex", "CurAge", "isDead")]))
})))
afftime_df$TimeAllCancers = apply(afftime_df[,grep("Time*", names(afftime_df))], 
                                  1, min)

save(mut_df, afftime_df, file="results/output/est_pen_dfs.rData")

###############################################################################

cancer_map = PanelPRO:::CANCER_NAME_MAP$short
names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long

# Initialize array for storing direct estimates
combs = expand.grid(ALLGENES = list(0, 1), Sex = c(0, 1))
est_pen = array(dim = c(3, 2, 2, 18, 94), 
                dimnames = list(
                  value = c("counts", "prob", "smoothed"), 
                  ALLGENES = c(0, 1), 
                  Sex = c(0, 1), 
                  Cancer = c(cancers, "AllCancers"), 
                  Age = 1:94))

# Iterate through all of the gene/sex combinations
for (i in 1:nrow(combs)) {
  # Number of cancer cases at each age interval
  counts = apply(
    afftime_df[paste0("Time", c(cancer_map[cancers], "AllCancers"))], 
    2, function(x) {
    table(factor(x[mut_df$ALLGENES==combs$ALLGENES[i] & 
                     afftime_df$Sex==combs$Sex[i]], levels = 1:95))
  })[-95,]
  
  # Divide by the number of people with the corresponding carrier status and sex
  prob = counts / (sum(mut_df$ALLGENES==combs$ALLGENES[i] & 
                         afftime_df$Sex==combs$Sex[i])) 
  
  # Apply smoothing spline
  smoothed = apply(prob, 2, function(x) {smooth.spline(x, spar = 0.6)$y})
  # Replace negative values with 0
  smoothed[smoothed < 0] = 0
  
  est_pen["counts", as.character(combs$ALLGENES[i]), as.character(combs$Sex[i]), ,] = t(counts)
  est_pen["prob", as.character(combs$ALLGENES[i]), as.character(combs$Sex[i]), ,] = t(prob)
  est_pen["smoothed", as.character(combs$ALLGENES[i]), as.character(combs$Sex[i]), ,] = t(smoothed)
}

save(est_pen, file="estimated_parameters.rData")
