library(tidyverse)
library(gtsummary)

# Load HCP and MGH families
load("../../Clinical_validation/hcp/data/pp21_other_cancers/hcp_families.rData")
load("../../Clinical_validation/mgh/data/pp21_other_cancers/mgh_families.rData")

# Only keep MGH families with 21 genes tested
panel_id_21 = which(sapply(panel_list, length) == 21)
mgh_families_PP = mgh_families_PP[proband_muts$PanelID == panel_id_21]

# Remove families that raised a model error
load("./hcp/results/diagnostics/diagnostics.rData")
hcp_families_PP = hcp_families_PP[!failed_families]
load("./mgh/results/diagnostics/diagnostics.rData")
mgh_families_PP = mgh_families_PP[!failed_families]

# Initialize list for storing families
all_dat = list(HCP = hcp_families_PP, MGH = mgh_families_PP)

# PanelPRO-21 genes, cancers, and tumor markers
genes = sort(c("ATM", "BARD1", "BRCA1", "BRCA2", "BRIP1", "CDH1", 
               "CDK4", "CDKN2A", "CHEK2", "EPCAM",  "MLH1", "MSH2", 
               "MSH6", "NBN", "PALB2", "PMS2", "PTEN", "RAD51C", 
               "RAD51D", "STK11", "TP53"))
cancers = c("BRA", "BC", "COL", "ENDO", "GAS", 
            "HEP", "KID", "LEUK", "MELA", "OC", "OST", 
            "PANC", "PROS", "SI", "STS", "THY", "UB")
markers = c("ER", "CK14", "CK5.6", 
            "PR", "HER2", "MSI")

# Create data frame with all variables, at the counselee level, 
# for both cohorts
my_df = do.call(rbind, lapply(names(all_dat), function(x) {
  data.frame(# Number of counselees
    #N = length(x), 
    # Cohort
    Cohort = x, 
    # Mean number of people in each family
    MeanSize = sapply(all_dat[[x]], nrow), 
    # Counselee sex
    Male = sapply(all_dat[[x]], function(fam) {
      fam$Sex[fam$isProband==1]
    }),  
    # Counselee age
    Age = sapply(all_dat[[x]], function(fam) {
      fam$CurAge[fam$isProband==1]
    }), 
    # Counselee ancestry
    Ancestry = sapply(all_dat[[x]], function(fam) {
      fam$Ancestry[fam$isProband==1]
    }), 
    # Counselee race
    Race = sapply(all_dat[[x]], function(fam) {
      fam$race[fam$isProband==1]
    }), 
    # Total number of counselees who are carriers
    AnyCarrier = sapply(all_dat[[x]], function(fam){
      any(as.numeric(fam[fam$isProband==1, genes]) == 1, na.rm = TRUE)
    }),
    # Total number of counselees who are affected for cancers in the model
    AnyCancer = sapply(all_dat[[x]], function(fam){
      any(fam[fam$isProband==1,paste0("isAff", cancers)] > 0)
    }),
    # Total number of relatives who are affected for cancers in the model
    AnyRelCancer = sapply(all_dat[[x]], function(fam){
      sum(apply(fam[,paste0("isAff", cancers)] > 0, 1, any))
    })
  )
}))

# Set default ancestry and race to missing values
my_df$Ancestry[is.na(my_df$Ancestry)] = "nonAJ"
my_df$Race[is.na(my_df$Race)] = "All_Races"

# Nicer factor levels for race
my_df = my_df %>% 
  mutate(
    Race = fct_recode(Race, "All Races" = "All_Races"))

# Create a "Table 1"
my_df %>% 
  tbl_summary(
    by = "Cohort", 
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 2,
      all_categorical() ~ c(0, 2)), 
    label = list(
      MeanSize = "Family size, mean (SD)", 
      Male = "Male counselees, n (%)", 
      Age = "Counselee age in years, mean (SD)", 
      Ancestry = "Counselee ancestry, n (%)", 
      Race = "Counselee race, n (%)", 
      AnyCarrier = "Counselee carriers of a pathogenic variant on any gene, n (%)", 
      AnyCancer = "Counselees diagnosed with any cancer, n (%)", 
      AnyRelCancer = "Family members diagnosed with any cancer, mean (SD)")) %>% 
  as_kable(format = "latex")
