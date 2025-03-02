library(abind)

K = 100
# Read in bootstrap diagnostics
all_diagnostics_boot = NULL
for (k in 1:K){
  load(paste0("results/diagnostics/diagnostics_boot", k, ".rData"))
  all_diagnostics_boot = c(all_diagnostics_boot, diagnostics_boot)
}

# Pull out and reshape relevant diagnostics
reshaped_boot_diagnostics = abind(all_diagnostics_boot, along = 3)

# Get summary statistics of bootstrap diagnostics
boot_diagnostics_summary = 
  abind(apply(reshaped_boot_diagnostics, c(1,2), function(y){summary(y)[1:6]}), 
        SD=array(apply(reshaped_boot_diagnostics, c(1,2), sd), dim=c(1,dim(reshaped_boot_diagnostics)[1:2])), 
        CI95lo=array(apply(reshaped_boot_diagnostics, c(1,2), quantile, probs=0.025, na.rm=TRUE), 
                     dim=c(1,dim(reshaped_boot_diagnostics)[1:2])), 
        CI95hi=array(apply(reshaped_boot_diagnostics, c(1,2), quantile, probs=0.975, na.rm=TRUE), 
                     dim=c(1,dim(reshaped_boot_diagnostics)[1:2])), 
        along=1)

# Save bootstraps
save(reshaped_boot_diagnostics, 
     file="results/diagnostics/boot_diagnostics.rData")

# Save summary statistics
save(boot_diagnostics_summary, 
     file="results/diagnostics/boot_diagnostics_summary.rData")
