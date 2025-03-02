# Discrimination (AUC)
library(pROC)

# Calibration (does it sum up to the number of carriers?)
# Ratio of the sum of observed to the sum of predicted 
# Drops NAs
# x = vector of observed values
# pred = vector of predicted values
calibration = function(x, pred){
  na_idx = is.na(x) | is.na(pred)
  return((sum(x[!na_idx])/sum(pred[!na_idx])))
}

# Mean squared error
# Drops NAs
# x = vector of observed values
# pred = vector of predicted values
mse = function(x, pred) {
  return(mean((x-pred)^2, na.rm=TRUE))
}

# Obtain all three diagnostics
# x = vector of observed values
# pred = vector of predicted values
# title = character string to include in the plot title besides the AUC
# ... additional arguments to pass to roc.plot
get_diagnostics = function(x, pred) {
  return(c(AUC = tryCatch(auc(roc(x, pred)), error=function(err) NA),
           "calib. O/E" = calibration(x, pred),
           MSE = mse(x, pred)))
}


get_all_diagnostics = function(labels, probs, return_probs = TRUE) {
  
  # Diagnostics 
  out = get_diagnostics(labels, probs)
  
  # Output with estimated probabilities
  if (return_probs==TRUE) {
    out = list(diagnostics=out, 
               probs=probs, 
               labels=labels)
  } 
  
  return(out)
}
