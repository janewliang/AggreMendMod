# Load libraries
library(PanelPRO)
library(abind)

# Add OtherCancers to the list of acceptable cancers
NEW_CANCER_TYPES = c(PanelPRO:::CANCER_TYPES, "OtherCancers")
NEW_CANCER_TYPES = NEW_CANCER_TYPES[!duplicated(NEW_CANCER_TYPES)]
assignInNamespace("CANCER_TYPES", NEW_CANCER_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_CANCER_NAME_MAP = list(short = c(PanelPRO:::CANCER_NAME_MAP$short, "OtherCancers"), 
                           long = c(PanelPRO:::CANCER_NAME_MAP$long, "OtherCancers"))
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
mult_fact = 3

# New database
new_db = PanelPRODatabase

# Set up OtherCancers penetrance for each gene
new_db$Penetrance = abind(new_db$Penetrance, 
                          OtherCancers=array(0, dim=dim(new_db$Penetrance[1,,,,,])), 
                          along=1)

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
                   OtherCancers=new_db$DOC["Brain",,,,], 
                   along=1)
names(dimnames(new_db$DOC)) = c("Cancer", "Gene", "Race", "Sex", "Age")

# Riskmod
new_db$Riskmod = abind(new_db$Riskmod, 
                       OtherCancers=array(1, dim=dim(new_db$Riskmod[1,,,,,])), 
                       along=1)
names(dimnames(new_db$Riskmod)) = c("Cancer", "Gene", "Intervention", 
                                    "Sex", "IntervAge", "DataType")
