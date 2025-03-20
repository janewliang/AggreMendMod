library(tidyverse)
library(reshape2)
library(abind)

# Set background to be white for all ggplots
theme_set(theme_bw() + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())) 

# Read in diagnostics
load("./results/diagnostics/diagnostics.rData")
load("./results/diagnostics/boot_diagnostics_summary.rData")

# Set up data for plotting
plot_df = list(melt(diagnostics, value.name = "est"), 
               melt(boot_diagnostics_summary["CI95lo",,], value.name = "lo"), 
               melt(boot_diagnostics_summary["CI95hi",,], value.name = "hi")) %>% 
  reduce(left_join) %>% 
  rename(metric = Var1, mod = Var2) %>% 
  mutate(metric = factor(fct_recode(metric, "O/E" = "calib. O/E"), 
                         levels = c("AUC", "O/E", "MSE")), 
         mod_type = factor(rep(rep(c("Fam3PRO", "Aggregate"), each = 3), 3), 
                           levels = c("Fam3PRO", "Aggregate")), 
         perc_misreport = rep(c(10, 30, 50), each = 6), 
         mod = factor(paste0(mod_type, " (", perc_misreport, "%)"), 
                      levels = c(paste0(c("Fam3PRO", "Aggregate"), " (", 
                                        rep(c(10, 30, 50), each = 2), "%)"))))

# Plot diagnostic metrics
plot_df %>% 
  ggplot(aes(x = mod, y = est, color = mod_type, shape = mod_type)) + 
  geom_pointrange(aes(ymin = lo, ymax = hi)) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) + 
  scale_shape_manual(values = c(17, 16)) + 
  facet_wrap(~metric, scales = "free", strip.position = "left") + 
  geom_vline(xintercept = c(2.5, 4.5), linetype = 2, color = "grey") + 
  xlab("") + ylab("") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"), 
        strip.placement = "outside", 
        panel.spacing = unit(0.5, "lines"), 
        plot.margin = unit(c(0.5, 0.5, -0.5, -1), "lines"), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("./results/diagnostics/diagnostics.png", 
       width = 6, height = 4, dpi = 300)

# Tables of diagnostic metrics
plot_df %>% 
  mutate(value = paste0(round(est, 3), 
                        " (", round(lo, 3), ", ", round(hi, 3), ")")) %>% 
  dplyr::select(metric, mod, value) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  knitr::kable(format = "latex")
