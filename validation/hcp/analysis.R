library(abind)

metrics = c("AUC"=1, "E/O"=2, "MSE"=5)
mods = -3
multiplying_factors = seq(1.2, 3, by = 0.3)

# Read in diagnostics
load("results/diagnostics/diagnostics.rData")

# Pull out and reshape relevant diagnostics
reshaped_diagnostics = do.call(abind, 
                               c(lapply(diagnostics, function(x) {
                                 x[metrics,mods]
                               }), along = 3))

load("results/diagnostics/boot_diagnostics_summary.rData")
# Put all bootstrap confidence intervals together
# Lower CI
reshaped_diagnostics_CI95lo = 
  lapply(metrics, 
         function(x){
           sapply(boot_diagnostics_summary, function(y){
             y["CI95lo",x,mods]
           })
         })
# Upper CI
reshaped_diagnostics_CI95hi = 
  lapply(metrics, 
         function(x){
           sapply(boot_diagnostics_summary, function(y){
             y["CI95hi",x,mods]
           })
         })


plot_dat = lapply(1:dim(reshaped_diagnostics)[3], function(i) {
  out = abind(reshaped_diagnostics[,,i], 
              t(sapply(reshaped_diagnostics_CI95lo, function(x) {x[,i]})), 
              t(sapply(reshaped_diagnostics_CI95hi, function(x) {x[,i]})), along=3)
  dimnames(out)[[2]] = c("PanelPRO-22\n(risk mod.)", 
                         "PanelPRO-22", 
                         "Aggregate\n(17 cancers)", 
                         rep("Aggregate\n(17+ cancers)", length(multiplying_factors)))
  dimnames(out)[[3]] = c("Est", "CI95lo", "CI95hi")
  return(out)
})

# Make plots comparing diagnostics
my_cols = c("#D55E00", "#D55E00", 
            rep("#0072B2", times = 1+ length(multiplying_factors)))
my_pch = c(c(15, 17, 16), rep(1, times = length(multiplying_factors)))

# y-axis limits
ymin = sapply(plot_dat, function(x){
  apply(x, 1, min, na.rm=TRUE)
})
ymax = sapply(plot_dat, function(x){
  apply(x, 1, max, na.rm=TRUE)
})

png("results/diagnostics/diagnostics_3_panel_%01d.png", width=700, height=430)
par(mfrow=c(1,3), mar=c(7.5, 5, 2.8, 1), xpd=TRUE)
invisible(sapply(1:length(plot_dat), function(i) {
  plot_dat_i = plot_dat[[i]][,c(1:3, 5),]
  sapply(dimnames(plot_dat_i)[[1]], function(met) {
    plot(plot_dat_i[met,,"Est"], 
         xlim=c(0.5,4.5),
         ylim=c(ymin[met,i], ymax[met,i]), xaxt="n", 
         main=met, xlab="", ylab=met, pch=my_pch, col=my_cols, 
         cex=2, cex.lab=1.8, cex.axis=1.8, cex.main=2)
    segments(x0=1:dim(plot_dat_i)[2], 
             y0=plot_dat_i[met,,"CI95lo"], 
             y1=plot_dat_i[met,,"CI95hi"], 
             col=my_cols)
    # x-axis labels
    axis(1, 1:dim(plot_dat_i)[2], FALSE)
    text(x=1:dim(plot_dat_i)[2], 
         y=par("usr")[3]-(ymax[met,i]-ymin[met,i])/20, 
         srt=45, adj=1, 
         labels=dimnames(plot_dat_i)[[2]], 
         xpd=TRUE, cex=1.4)
    
    # Reference line for AUC and calibration
    if(met=="AUC" || met=="E/O" || met=="O/E") {
      segments(x0=0.65, x1=dim(plot_dat_i)[2]+0.35,
               y0=1, col="grey", lty=2)
    } else if (met=="E-O") {
      segments(x0=0.65, x1=dim(plot_dat_i)[2]+0.34,
               y0=0, col="grey", lty=2)
    }
  })
}))
dev.off()

png("results/diagnostics/RMSE.png", width=600, height=430)
par(mar=c(5, 5, 2.8, 1), xpd=TRUE)
rmse_dat = sqrt(plot_dat[[1]]["MSE",-(1:3),])
plot(rmse_dat[,"Est"], xaxt="n", 
     main="RMSE", xlab="Multiplying factor", ylab="RMSE", 
     pch=my_pch[-(1:3)], col=my_cols[-(1:3)], 
     cex=2, cex.lab=1.8, cex.axis=1.8, cex.main=2)
# x-axis labels
axis(1, 1:length(rmse_dat[,"Est"]), multiplying_factors, cex.axis=1.4)
dev.off()


# Tables of diagnostic metrics
lapply(names(metrics), function(x) {
  tab = data.frame(c("PanelPRO-22 (risk mod.)", 
                     "PanelPRO-22", 
                     "Aggregate (17 cancers)", 
                     rep("Aggregate (17+ cancers)", 
                         length(multiplying_factors))), 
                   c(rep("", 3), multiplying_factors), 
                   plot_dat[[1]][x,,])
  names(tab) = c("Model", "Mult. fact.", "Estimate", 
                 "Bootstrap 2.5%", "Bootstrap 97.5%")
  rownames(tab) = NULL
  knitr::kable(tab, digits = 5, format = "latex")
})

