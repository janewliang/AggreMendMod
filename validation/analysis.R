library(tidyverse)
library(reshape2)
library(abind)

# Set background to be white for all ggplots
theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())) 

# Read in diagnostics
# HCP
load("hcp/results/diagnostics/diagnostics.rData")
load("hcp/results/diagnostics/boot_diagnostics_summary.rData")
diagnostics_hcp = diagnostics
boot_diagnostics_summary_hcp = boot_diagnostics_summary
# MGH
load("mgh/results/diagnostics/diagnostics.rData")
load("mgh/results/diagnostics/boot_diagnostics_summary.rData")
diagnostics_mgh = diagnostics
boot_diagnostics_summary_mgh = boot_diagnostics_summary

# Set up data for plotting
# HCP
plot_df_hcp = 
  list(melt(diagnostics_hcp, value.name = "est"), 
       melt(boot_diagnostics_summary_hcp["CI95lo",,], value.name = "lo"), 
       melt(boot_diagnostics_summary_hcp["CI95hi",,], value.name = "hi")) %>% 
  reduce(left_join) %>% 
  rename(metric = Var1, mod = Var2) %>% 
  mutate(metric = factor(fct_recode(metric, "O/E" = "calib. O/E"), 
                         levels = c("AUC", "O/E", "MSE")), 
         mod = factor(fct_recode(mod, "Fam3PRO" = "PP21", "Aggregate" = "Agg"), 
                      levels = c("Fam3PRO", "Aggregate")))
# MGH
plot_df_mgh = 
  list(melt(diagnostics_mgh, value.name = "est"), 
       melt(boot_diagnostics_summary_mgh["CI95lo",,], value.name = "lo"), 
       melt(boot_diagnostics_summary_mgh["CI95hi",,], value.name = "hi")) %>% 
  reduce(left_join) %>% 
  rename(metric = Var1, mod = Var2) %>% 
  mutate(metric = factor(fct_recode(metric, "O/E" = "calib. O/E"), 
                         levels = c("AUC", "O/E", "MSE")), 
         mod = factor(fct_recode(mod, "Fam3PRO" = "PP21", "Aggregate" = "Agg"), 
                      levels = c("Fam3PRO", "Aggregate")))
# Combined
plot_df = rbind(
  data.frame(data = "HCP", plot_df_hcp), 
  data.frame(data = "MGH", plot_df_mgh)
) %>% 
  mutate(
    data_mod = factor(paste(data, "-", mod), 
                      levels = c("HCP - Fam3PRO", 
                                 "HCP - Aggregate", 
                                 "MGH - Fam3PRO", 
                                 "MGH - Aggregate")))

# Plot diagnostic metrics
plot_df %>% 
  ggplot(aes(x = data_mod, y = est, color = mod, shape = mod)) + 
  geom_pointrange(aes(ymin = lo, ymax = hi)) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  scale_shape_manual(values = c(17, 16)) + 
  facet_wrap(~metric, scales = "free", strip.position = "left") + 
  geom_vline(xintercept = 2.5, linetype = 2, color = "grey") + 
  xlab("") + ylab("") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"), 
        strip.placement = "outside", 
        panel.spacing = unit(0.5, "lines"), 
        plot.margin = unit(c(0.5, 0.5, -0.5, -1), "lines"), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("diagnostics.png", 
       width = 6, height = 4, dpi = 300)

# Tables of diagnostic metrics
plot_df %>% 
  mutate(value = paste0(round(est, 3), 
                        " (", round(lo, 3), ", ", round(hi, 3), ")")) %>% 
  dplyr::select(Data = data, metric, Model = mod, value) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  knitr::kable(format = "latex")
