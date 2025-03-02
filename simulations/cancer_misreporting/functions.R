# Load libraries
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

# Add AllCancers to the list of acceptable cancers
NEW_CANCER_TYPES = c(PanelPRO:::CANCER_TYPES, "AllCancers")
NEW_CANCER_TYPES = NEW_CANCER_TYPES[!duplicated(NEW_CANCER_TYPES)]
assignInNamespace("CANCER_TYPES", NEW_CANCER_TYPES, 
                  ns="PanelPRO", pos="package:PanelPRO")
NEW_CANCER_NAME_MAP = list(short = c(PanelPRO:::CANCER_NAME_MAP$short, "AllCancers"), 
                           long = c(PanelPRO:::CANCER_NAME_MAP$long, "AllCancers"))
NEW_CANCER_NAME_MAP = list(short = NEW_CANCER_NAME_MAP$short[!duplicated(NEW_CANCER_NAME_MAP$short)], 
                           long = NEW_CANCER_NAME_MAP$long[!duplicated(NEW_CANCER_NAME_MAP$long)])
assignInNamespace("CANCER_NAME_MAP", NEW_CANCER_NAME_MAP, 
                  ns="PanelPRO", pos="package:PanelPRO")

###############################################################################

# Load estimated penetrances
load("../../est_pen/pp21/estimated_parameters.rData")

# Estimate ALLGENES allele frequency
af_idx = PanelPRO:::formatGeneNames(rownames(PanelPRODatabase$AlleleFrequency), 
                                    format = "only_gene") %in% genes
af = sum(PanelPRODatabase$AlleleFrequency[af_idx,"nonAJ"]) + 
  sum(apply(combn(PanelPRODatabase$AlleleFrequency[af_idx,"nonAJ"], 2), 2, prod))
# 0.01137739

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

# Set up AllCancers penetrance for each gene
new_db$Penetrance = abind(new_db$Penetrance, 
                          AllCancers=array(0, dim=dim(new_db$Penetrance[1,,,,,])), 
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

names(dimnames(new_db$Penetrance)) = c("Cancer", "Gene", "Race", 
                                       "Sex", "Age", "PenetType")

# DOC
new_db$DOC = abind(new_db$DOC, 
                   ALLGENES=new_db$DOC[,"SEER",,,], 
                   along=2)
new_db$DOC = abind(new_db$DOC, 
                   AllCancers=new_db$DOC["Brain",,,,], 
                   along=1)
names(dimnames(new_db$DOC)) = c("Cancer", "Gene", "Race", "Sex", "Age")

# Riskmod
new_db$Riskmod = abind(new_db$Riskmod, 
                       ALLGENES=array(1, dim=dim(new_db$Riskmod[,1,,,,])), 
                       along=2)
new_db$Riskmod = abind(new_db$Riskmod, 
                       AllCancers=array(1, dim=dim(new_db$Riskmod[1,,,,,])), 
                       along=1)
names(dimnames(new_db$Riskmod)) = c("Cancer", "Gene", "Intervention", 
                                    "Sex", "IntervAge", "DataType")

# Germline testing
new_db$GermlineTesting = abind(new_db$GermlineTesting, 
                               ALLGENES=rep(1, dim(new_db$GermlineTesting)[2]), 
                               along=1)

################################################################################

# PanelPRO-21 
PanelPRO21 = function(...) {
  PanelPRO::PanelPRO(
    genes = c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDKN2A", 
              "CDK4", "CHEK2", "EPCAM", "MLH1", "MSH2", "MSH6", "NBN", 
              "PALB2", "PMS2", "PTEN", "RAD51C", "RAD51D", "STK11", "TP53"), 
    cancers = c("Brain", "Breast", "Colorectal", "Endometrial", "Gastric", 
                "Kidney", "Leukemia", "Melanoma", "Ovarian", "Osteosarcoma", 
                "Pancreas", "Prostate", "Small Intestine", 
                "Soft Tissue Sarcoma", "Thyroid", "Urinary Bladder", 
                "Hepatobiliary", "Contralateral"), 
    ...)
}

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

# Aggregate
run_agg = function(...) {
  PanelPRO::PanelPRO(genes="ALLGENES", 
                     cancers="AllCancers", 
                     ...)
}

# Function to re-assign cancer types and ages to a random subset of people
swapCancerTypes = function(fam, prob) {
  
  # List of cancers
  female_cancers = c("BRA", "BC", "COL", "ENDO", "GAS", # No prostate
                     "KID", "LEUK", "MELA", "OC", "OST", 
                     "PANC", "SI", "STS", "THY", "UB", "HEP")
  male_cancers = c("BRA", "BC", "COL", "GAS", "KID", # No endometrial, ovarian
                   "LEUK", "MELA", "OST", "PANC", "PROS", 
                   "SI", "STS", "THY", "UB", "HEP")
  
  # Identify individuals who are affected for cancer
  cancer_idx = which(fam$isAffAny == 1)
  
  # Randomly sample prob% to re-assign cancers
  if (length(cancer_idx) > 0) {
    jumble_cancer_idx = cancer_idx[rbinom(length(cancer_idx), 1, prob)==1]
    
    # If there are people who should have missing cancers, remove them
    if (length(jumble_cancer_idx) > 0) {
      for (idx in jumble_cancer_idx) {
        if (fam$Sex[idx] == 0) { # Females
          # Indices to re-order cancers
          new_cancer_order = sample(length(female_cancers))
          
          # Re-assign which cancers that person has
          fam[idx, paste0("isAff", female_cancers)] = 
            fam[idx, paste0("isAff", female_cancers)][new_cancer_order]
          fam[idx, paste0("Age", female_cancers)] = 
            fam[idx, paste0("Age", female_cancers)][new_cancer_order]
          fam[idx, paste0("Time", female_cancers)] = 
            fam[idx, paste0("Time", female_cancers)][new_cancer_order]
        } else { # Males
          # Indices to re-order cancers
          new_cancer_order = sample(length(male_cancers))
          
          # Re-assign which cancers that person has
          fam[idx, paste0("isAff", male_cancers)] = 
            fam[idx, paste0("isAff", male_cancers)][new_cancer_order]
          fam[idx, paste0("Age", male_cancers)] = 
            fam[idx, paste0("Age", male_cancers)][new_cancer_order]
          fam[idx, paste0("Time", male_cancers)] = 
            fam[idx, paste0("Time", male_cancers)][new_cancer_order]
        }
      }
    } 
  }
  
  return(fam)
}

# Function to obtain posterior probabilities
run_models = function(fam_PP) {
  
  # PanelPRO-21
  probs_pp21 = extract_probs(PanelPRO21(fam_PP, database = new_db,
                                        parallel = FALSE, 
                                        net = TRUE, max.mut = 2, 
                                        allow.intervention = FALSE))
  
  # Aggregate
  probs_agg = extract_probs(run_agg(fam_PP, database = new_db,
                                    parallel = FALSE, 
                                    net = TRUE, 
                                    allow.intervention = FALSE))
  
  # Carrier probabilities
  out = rbind(c(probs_pp21[1], 1-probs_pp21[1]), 
              probs_agg)
  rownames(out) = c("pp21", "agg")
  colnames(out) = c("noncarrier", "ALLGENES")
  
  return(out)
}
