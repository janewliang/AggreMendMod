# Load libraries
library(PanelPRO)
library(abind)

# Add ALLGENES gene to the list of acceptable genes
NEW_GENE_TYPES = c(PanelPRO:::GENE_TYPES, "ALLGENES")
NEW_GENE_TYPES = NEW_GENE_TYPES[!duplicated(NEW_GENE_TYPES)]
assignInNamespace("GENE_TYPES", NEW_GENE_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_ALL_GENE_VARIANT_TYPES = c(PanelPRO:::ALL_GENE_VARIANT_TYPES, 
                               ALLGENES = "ALLGENES")
NEW_ALL_GENE_VARIANT_TYPES = NEW_ALL_GENE_VARIANT_TYPES[!duplicated(NEW_ALL_GENE_VARIANT_TYPES)]
assignInNamespace("ALL_GENE_VARIANT_TYPES", NEW_ALL_GENE_VARIANT_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_DEFAULT_VARIANTS = c(PanelPRO:::DEFAULT_VARIANTS, 
                         ALLGENES = "ALLGENES")
NEW_DEFAULT_VARIANTS = NEW_DEFAULT_VARIANTS[!duplicated(NEW_DEFAULT_VARIANTS)]
assignInNamespace("DEFAULT_VARIANTS", NEW_DEFAULT_VARIANTS, 
                  ns="PanelPRO", pos="package:PanelPRO")

# Add AllCancers and OtherCancers to the list of acceptable cancers
NEW_CANCER_TYPES = c(PanelPRO:::CANCER_TYPES, "AllCancers", "OtherCancers")
NEW_CANCER_TYPES = NEW_CANCER_TYPES[!duplicated(NEW_CANCER_TYPES)]
assignInNamespace("CANCER_TYPES", NEW_CANCER_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_CANCER_NAME_MAP = list(short = c(PanelPRO:::CANCER_NAME_MAP$short, "AllCancers", "OtherCancers"), 
                           long = c(PanelPRO:::CANCER_NAME_MAP$long, "AllCancers", "OtherCancers"))
NEW_CANCER_NAME_MAP = list(short = NEW_CANCER_NAME_MAP$short[!duplicated(NEW_CANCER_NAME_MAP$short)], 
                           long = NEW_CANCER_NAME_MAP$long[!duplicated(NEW_CANCER_NAME_MAP$long)])
assignInNamespace("CANCER_NAME_MAP", NEW_CANCER_NAME_MAP, 
                  ns="PanelPRO", pos="package:PanelPRO")

###############################################################################

# 17 cancers to consider
cancers = c("Brain", "Breast", "Colorectal", "Endometrial", 
            "Gastric", "Hepatobiliary", "Kidney", "Leukemia", 
            "Melanoma", "Osteosarcoma", "Ovarian", "Pancreatic", 
            "Prostate", "SmallIntestine", "SoftTissueSarcoma", 
            "Thyroid", "UrinaryBladder")
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
# Raw SEER counts for each age interval
SEER_f_counts = sapply(c("All", female_cancers), function(canc){
  as.numeric(gsub(",", "", 
                  read.table(paste0("../../SEER/SEER_1517_Female_", 
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

# Male
# Raw SEER counts for each age interval
SEER_m_counts = sapply(c("All", male_cancers), function(canc){
  as.numeric(gsub(",", "", 
                  read.table(paste0("../../SEER/SEER_1517_Male_", 
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

###############################################################################

# Multiplying factor for OtherCancers carrier penetrance
mult_fact = 1.5

# Load estimated penetrances
load(paste0("../../est_pen/oc", mult_fact, "/estimated_parameters.rData"))

af = sum(PanelPRODatabase$AlleleFrequency[,"nonAJ"]) + 
  sum(apply(combn(PanelPRODatabase$AlleleFrequency[,"nonAJ"], 2), 2, prod))
# 0.01145489

# New database
new_db = PanelPRODatabase

## Hypothetical ALLGENES
# Add ALLGENES allele frequency
new_db$AlleleFrequency = abind(new_db$AlleleFrequency, 
                               ALLGENES=rep(af, dim(new_db$AlleleFrequency)[2]), 
                               along=1)
names(dimnames(new_db$AlleleFrequency)) = c("Gene", "Ancestry")

# Set up ALLGENES penetrance for each cancer
new_db$Penetrance = abind(new_db$Penetrance, 
                          ALLGENES=array(0, dim=dim(new_db$Penetrance[,1,,,,])), 
                          along=2)

# Set up AllCancers and OtherCancers penetrance for each gene
new_db$Penetrance = abind(new_db$Penetrance, 
                          AllCancers=array(0, dim=dim(new_db$Penetrance[1,,,,,])), 
                          OtherCancers=array(0, dim=dim(new_db$Penetrance[1,,,,,])), 
                          along=1)

# Fill in AllCancers penetrances estimated based on simulations for each cancer
# SEER penetrances
new_db$Penetrance["AllCancers","SEER","All_Races","Female",,"Net"] = 
  est_pen["smoothed", "0", "0", "AllCancers", ] * (1-af)^2 + 
  est_pen["smoothed", "1", "0", "AllCancers", ] * 2*af*(1-af)
new_db$Penetrance["AllCancers","SEER","All_Races","Male",,"Net"] = 
  est_pen["smoothed", "0", "1", "AllCancers", ] * (1-af)^2 + 
  est_pen["smoothed", "1", "1", "AllCancers", ] * 2*af*(1-af)
# ALLGENES penetrances
new_db$Penetrance["AllCancers","ALLGENES","All_Races","Female",,"Net"] = 
  est_pen["smoothed", "1", "0", "AllCancers", ]
new_db$Penetrance["AllCancers","ALLGENES","All_Races","Male",,"Net"] =
  est_pen["smoothed", "1", "1", "AllCancers", ]

# Fill in OtherCancers penetrances with SEER/carrier penetrances generated by 
# multiplying the mult_fact
for (gene in dimnames(new_db$Penetrance)[[2]]) {
  if (gene == "SEER") { # SEER/population
    # SEER estimate for OtherCancers
    new_db$Penetrance["OtherCancers",gene,"All_Races","Female",,"Net"] = 
      OtherCancers_f
    new_db$Penetrance["OtherCancers",gene,"All_Races","Male",,"Net"] = 
      OtherCancers_m
  } else { # Genes/carriers
    # SEER estimate for OtherCancers, multiplied by mult_fact
    new_db$Penetrance["OtherCancers",gene,"All_Races","Female",,"Net"] = 
      OtherCancers_f * mult_fact
    new_db$Penetrance["OtherCancers",gene,"All_Races","Male",,"Net"] = 
      OtherCancers_m * mult_fact
  }
}

names(dimnames(new_db$Penetrance)) = c("Cancer", "Gene", "Race", 
                                       "Sex", "Age", "PenetType")

# DOC
new_db$DOC = abind(new_db$DOC, 
                   ALLGENES=new_db$DOC[,"SEER",,,], 
                   along=2)
new_db$DOC = abind(new_db$DOC, 
                   AllCancers=new_db$DOC["Brain",,,,], 
                   OtherCancers=new_db$DOC["Brain",,,,], 
                   along=1)
names(dimnames(new_db$DOC)) = c("Cancer", "Gene", "Race", "Sex", "Age")

# Riskmod
new_db$Riskmod = abind(new_db$Riskmod, 
                       ALLGENES=array(1, dim=dim(new_db$Riskmod[,1,,,,])), 
                       along=2)
new_db$Riskmod = abind(new_db$Riskmod, 
                       AllCancers=array(1, dim=dim(new_db$Riskmod[1,,,,,])), 
                       OtherCancers=array(1, dim=dim(new_db$Riskmod[1,,,,,])), 
                       along=1)
names(dimnames(new_db$Riskmod)) = c("Cancer", "Gene", "Intervention", 
                                    "Sex", "IntervAge", "DataType")

# Germline testing
new_db$GermlineTesting = abind(new_db$GermlineTesting, 
                               ALLGENES=rep(1, dim(new_db$GermlineTesting)[2]), 
                               along=1)

###############################################################################

# Load estimated penetrances and reweighing factors
load("../../est_pen/pp22/estimated_parameters.rData")

# New database
new_db_pp22 = new_db

# Fill in AllCancers penetrances estimated based on simulations for each cancer
# SEER penetrances
new_db_pp22$Penetrance["AllCancers","SEER","All_Races","Female",,"Net"] = 
  est_pen["smoothed", "0", "0", "AllCancers", ] * (1-af)^2 + 
  est_pen["smoothed", "1", "0", "AllCancers", ] * 2*af*(1-af)
new_db_pp22$Penetrance["AllCancers","SEER","All_Races","Male",,"Net"] = 
  est_pen["smoothed", "0", "1", "AllCancers", ] * (1-af)^2 + 
  est_pen["smoothed", "1", "1", "AllCancers", ] * 2*af*(1-af)
# ALLGENES penetrances
new_db_pp22$Penetrance["AllCancers","ALLGENES","All_Races","Female",,"Net"] = 
  est_pen["smoothed", "1", "0", "AllCancers", ]
new_db_pp22$Penetrance["AllCancers","ALLGENES","All_Races","Male",,"Net"] =
  est_pen["smoothed", "1", "1", "AllCancers", ]

names(dimnames(new_db_pp22$Penetrance)) = c("Cancer", "Gene", "Race", 
                                           "Sex", "Age", "PenetType")

################################################################################

# Extract posterior probabilities estimates as a named vector
# (ignores lower and upper CI)
extract_probs = function(res) {
  # Pull out posterior probabilities
  probs_df = res$posterior.prob[[1]]
  
  # Extract estimates only as a vector and name them using genotypes
  probs_vec = probs_df$estimate
  names(probs_vec) = probs_df$genes
  return(probs_vec)
}

run_agg = function(...) {
  PanelPRO::PanelPRO(genes="ALLGENES", 
                     cancers="AllCancers", 
                     ...)
}

# Function that modifies pedigree so that only first cancers are kept
firstCancerFam = function(fam) {
  
  # Get short cancer names
  cancer_map = PanelPRO:::CANCER_NAME_MAP$short
  names(cancer_map) = PanelPRO:::CANCER_NAME_MAP$long
  short_cancers = cancer_map[PanelPRO:::MODELPARAMS$PanPRO22$CANCERS]
  
  # Affection and age columns
  isaff_cols = paste0("isAff", short_cancers)
  age_cols = paste0("Age", short_cancers)
  
  # Iterate through all family members
  for (i in 1:nrow(fam)) {
    # Individual's cancer ages (for affected cancers)
    cancer_ages = fam[i,age_cols][fam[i,isaff_cols] == 1]
    
    if (length(cancer_ages) > 1) {
      # If there is more than one cancer age, attempt to pick the earliest one
      first_age_idx = which((cancer_ages) == min((cancer_ages), na.rm = T))
      
      if (length(first_age_idx) == 0) { # All cancer ages missing
        # Pick a random cancer to be first
        first_age_idx = sample(length(cancer_ages), 1)
      }  else if (length(first_age_idx) > 1) { # Multiple cancers at first age
        # Pick a random cancer to be first
        first_age_idx = sample(first_age_idx, 1)
      }
      
      # Set affection status for all other cancers to 0
      fam[i,isaff_cols][fam[i,isaff_cols] == 1][-first_age_idx] = 0
    }
  }
  return(fam)
}

# Function to obtain posterior probabilities
run_models = function(fam_PP) {
  
  # Only keep first cancers
  fam_PP_fc = firstCancerFam(fam_PP)
  
  # PanelPRO-22 (all cancers)
  probs_pp22_ac = extract_probs(PanelPRO::PanelPRO22(fam_PP, database = new_db,
                                                            parallel = FALSE, 
                                                            net = TRUE, max.mut = 2, 
                                                            allow.intervention = FALSE))
  
  # PanelPRO-22 (first cancers)
  probs_pp22_fc = extract_probs(PanelPRO::PanelPRO22(fam_PP_fc, database = new_db,
                                                            parallel = FALSE, 
                                                            net = TRUE, max.mut = 2, 
                                                            allow.intervention = FALSE))
  
  # Aggregate all cancers (17 cancers)
  probs_agg_noc = extract_probs(run_agg(fam_PP, database = new_db_pp22,
                                           parallel = FALSE, 
                                           net = TRUE, 
                                           allow.intervention = FALSE))
  
  # Aggregate all cancers (17 cancers plus "other cancers")
  probs_agg = extract_probs(run_agg(fam_PP, database = new_db,
                                            parallel = FALSE, 
                                           net = TRUE, 
                                           allow.intervention = FALSE))
  
  # Carrier probabilities
  out = rbind(c(probs_pp22_ac[1], 1-probs_pp22_ac[1]), 
                    c(probs_pp22_fc[1], 1-probs_pp22_fc[1]), 
                    probs_agg_noc,  probs_agg)
  rownames(out) = c("pp22_ac", "pp22_fc", "agg_noc", "agg")
  colnames(out) = c("noncarrier", "ALLGENES")
  
  return(out)
}
