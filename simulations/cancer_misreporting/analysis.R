library(abind)

metrics = c("AUC"=1, "E/O"=2, "MSE"=5)
mods = c(1, 3, 4, 6, 7, 9)

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
  dimnames(out)[[2]] = paste(rep(c("PanelPRO-22", 
                                   "Aggregate (17 cancers)"), 3), 
                             rep(c("(10%)", "(30%)", "(50%)"), each=2))
  dimnames(out)[[3]] = c("Est", "CI95lo", "CI95hi")
  return(out)
})

# Make plots comparing diagnostics
my_cols = rep(c("#D55E00", "#0072B2"), 3)
my_pch = rep(c(17, 16), 3)

# y-axis limits
ymin = sapply(plot_dat, function(x){
  apply(x, 1, min, na.rm=TRUE)
})
ymax = sapply(plot_dat, function(x){
  apply(x, 1, max, na.rm=TRUE)
})

png("results/diagnostics/diagnostics_3_panel_%01d.png", width=700, height=510)
par(mfrow=c(1,3), mar=c(16, 5, 2.8, 1), xpd=TRUE)
invisible(sapply(1:length(plot_dat), function(i) {
  plot_dat_i = plot_dat[[i]]
  sapply(dimnames(plot_dat_i)[[1]], function(met) {
    plot(plot_dat_i[met,,"Est"], 
         xlim=c(0.5,6.5),
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
         srt=55, adj=1, 
         labels=dimnames(plot_dat_i)[[2]], 
         xpd=TRUE, cex=1.4)
    
    # Reference lines separating % misreported
    segments(x0 = 2.5, y0 = par("usr")[3], y1 = par("usr")[4], 
             col="grey", lty=2)
    segments(x0 = 4.5, y0 = par("usr")[3], y1 = par("usr")[4], 
             col="grey", lty=2)
  })
}))
dev.off()


# Tables of diagnostic metrics
lapply(1:3, function(i) {
  lapply(names(metrics), function(x) {
    tab = data.frame(rep(c("PanelPRO-22", 
                           "Aggregate (17 cancers)"), 3), 
                     rep(c("10%", "30%", "50%"), each = 2), 
                     plot_dat[[i]][x,,])
    names(tab) = c("Model", "% Misreported", "Estimate", 
                   "Bootstrap 2.5%", "Bootstrap 97.5%")
    rownames(tab) = NULL
    return(knitr::kable(tab, format="latex", digits=5))
  })
})
