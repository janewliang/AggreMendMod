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


###############################################################################

# 17 cancers to consider
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", 
            "Gastric", "Hepatobiliary", "Kidney", "Leukemia", 
            "Melanoma", "Osteosarcoma", "Ovarian", "Pancreatic", 
            "Prostate", "SmallIntestine", "SoftTissueSarcoma", 
            "Thyroid", "UrinaryBladder")
# PanelPRO names
cancers_pp21 = c("Brain", "Breast", "Colorectal", "Endometrial", "Gastric", 
                 "Kidney", "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", 
                 "Pancreas", "Prostate", "Small Intestine", "Soft Tissue Sarcoma", 
                 "Thyroid", "Urinary Bladder", "Hepatobiliary")
female_cancers = c("Brain", "Breast", "Colorectal", # No prostate
                   "Endometrial", "Gastric", "Hepatobiliary", "Kidney", 
                   "Leukemia", "Melanoma", "Osteosarcoma", "Ovarian", 
                   "Pancreatic", "SmallIntestine", "SoftTissueSarcoma", 
                   "Thyroid", "UrinaryBladder")
male_cancers = c("Brain", "Breast", "Colorectal", # No endometrial, ovarian
                 "Gastric", "Hepatobiliary", "Kidney", "Leukemia", 
                 "Melanoma", "Osteosarcoma", "Pancreatic", "Prostate", 
                 "SmallIntestine", "SoftTissueSarcoma", "Thyroid", 
                 "UrinaryBladder")

# Female
# SEER net cancer penetrances
SEER_f = sapply(c("All", cancers), function(canc){
  if (canc %in% c("All", female_cancers)) {
    read.table(paste0("../SEER/SEER_1517_Female_", 
                      canc, "_Intermediate.txt"), nrows = 95)$V8[-1]
  } else {
    rep(0, 94)
  }
})
colnames(SEER_f) = c("AllCancers", sort(cancers_pp21))

# Raw SEER counts for each age interval
SEER_f_counts = sapply(c("All", female_cancers), function(canc){
  as.numeric(gsub(",", "", 
                  read.table(paste0("../SEER/SEER_1517_Female_", 
                                    canc, "_LT.txt"), nrows = 95)$V9[-1]))
})

# Calculate "Other Cancers" penetrance using raw counts 
# Need to re-scale each individual cancer's counts at each interval based on 
# 1 - cumulative sum of all other cancers in previous intervals 
# Idea is that the the counts for each interval are the people (out of 10000000)
# who developed the given cancer for the first time, but who may have gotten 
# other cancers in the previous intervals. 
# So assuming that the proportion of people who got cancer A in the interval 
# are equally likely as people who didn't get cancer to have developed the 
# other cancers, you can just subtract off [count cancer A] * proportion of 
# people who got cancer in previous intervals
SEER_f_counts_adj = SEER_f_counts[,-1]
OtherCancers_f = c(SEER_f_counts[1,"All"] - sum(SEER_f_counts[1,-1]), 
                   rep(0, nrow(SEER_f_counts)-1))
for (i in 2:length(OtherCancers_f)) {
  SEER_f_counts_adj[i,] = sapply(colnames(SEER_f_counts_adj), function(canc) {
    SEER_f_counts[i,canc] * 
      (1 - sum(SEER_f_counts_adj[1:(i-1),setdiff(colnames(SEER_f_counts_adj), canc)], 
               OtherCancers_f[1:(i-1)]) / 10000000)
  })
  OtherCancers_f[i] = SEER_f_counts[i,"All"] - sum(SEER_f_counts_adj[i,])
}
OtherCancers_f = OtherCancers_f / 10000000
SEER_f = cbind(SEER_f, OtherCancers = OtherCancers_f)


# Male
# SEER net cancer penetrances
SEER_m = sapply(c("All", cancers), function(canc){
  if (canc %in% c("All", male_cancers)) {
    read.table(paste0("../SEER/SEER_1517_Male_", 
                      canc, "_Intermediate.txt"), nrows = 95)$V8[-1]
  } else {
    rep(0, 94)
  }
})
colnames(SEER_m) = c("AllCancers", sort(cancers_pp21))

# Raw SEER counts for each age interval
SEER_m_counts = sapply(c("All", male_cancers), function(canc){
  as.numeric(gsub(",", "", 
                  read.table(paste0("../SEER/SEER_1517_Male_", 
                                    canc, "_LT.txt"), nrows = 95)$V9[-1]))
})

# Calculate "Other Cancers" penetrance using raw counts 
# Need to re-scale each individual cancer's counts at each interval based on 
# 1 - cumulative sum of all other cancers in previous intervals 
# Idea is that the the counts for each interval are the people (out of 10000000)
# who developed the given cancer for the first time, but who may have gotten 
# other cancers in the previous intervals. 
# So assuming that the proportion of people who got cancer A in the interval 
# are equally likely as people who didn't get cancer to have developed the 
# other cancers, you can just subtract off [count cancer A] * proportion of 
# people who got cancer in previous intervals
SEER_m_counts_adj = SEER_m_counts[,-1]
OtherCancers_m = c(SEER_m_counts[1,"All"] - sum(SEER_m_counts[1,-1]), 
                   rep(0, nrow(SEER_m_counts)-1))
for (i in 2:length(OtherCancers_m)) {
  SEER_m_counts_adj[i,] = sapply(colnames(SEER_m_counts_adj), function(canc) {
    SEER_m_counts[i,canc] * 
      (1 - sum(SEER_m_counts_adj[1:(i-1),setdiff(colnames(SEER_m_counts_adj), canc)], 
               OtherCancers_m[1:(i-1)]) / 10000000)
  })
  OtherCancers_m[i] = SEER_m_counts[i,"All"] - sum(SEER_m_counts_adj[i,])
}
OtherCancers_m = OtherCancers_m / 10000000
SEER_m = cbind(SEER_m, OtherCancers = OtherCancers_m)



# List of SEER penetrance data frames
SEER_dfs = list("0" = SEER_f, 
                "1" = SEER_m)


# Create a data frame with the multiplying factor between carrier and SEER
# penetrance for each cancer, at age each age
mult_fact_df = data.frame(
  do.call(cbind, 
          # Each sex
          lapply(c("0", "1"), function(sex) {
            # Each cancer
            sapply(dimnames(est_pen)$Cancer, function(cancer) {
              # Multiplying factor between carrier and SEER penetrance
              out = est_pen["smoothed", "1", sex, cancer,] / 
                SEER_dfs[[sex]][,cancer]
              # If either carrier/SEER penetrance was close to 0, set to 0
              out[est_pen["smoothed", "1", sex, cancer,] < 1e-05] = 0
              out[SEER_dfs[[sex]][,cancer] < 1e-05] = 0
              # Smooth the factors
              out = smooth.spline(out, spar = 0.6)$y
              # Replace negative values with 0
              out[out < 0] = 0
              out
            })
          })))
colnames(mult_fact_df) = paste(rep(c("Female", "Male"), 
                                   each = length(dimnames(est_pen)$Cancer)), 
                               dimnames(est_pen)$Cancer)
# Remove cancers that only apply to one gender
mult_fact_df[,colSums(mult_fact_df) == 0] = NULL

# For now, consider 3 scenerios
# Use female "All Cancers" for both genders
# Use male "All Cancers" for both genders
# Use female/male "All Cancers" for corresponding gender
mult_fact_all_df = list(
  Female = data.frame(Female = mult_fact_df[["Female AllCancers"]], 
                      Male = mult_fact_df[["Female AllCancers"]]), 
  Male = data.frame(Female = mult_fact_df[["Male AllCancers"]], 
                    Male = mult_fact_df[["Male AllCancers"]]), 
  Gendered = data.frame(Female = mult_fact_df[["Female AllCancers"]], 
                        Male = mult_fact_df[["Male AllCancers"]])
  )

# Save with estimated penetrances
save(est_pen, mult_fact_df, mult_fact_all_df, file="estimated_parameters.rData")
