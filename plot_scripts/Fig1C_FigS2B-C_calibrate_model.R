library(ggplot2)
library(ggpubr)

# FIGURE 1C
# Plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD.

file <- read.delim(file = '../output_files/calibration/gnomad_loci_obs_vs_scores.txt', header = TRUE, sep = "\t")

ggplot(data = file, aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 2) +
  labs(x = "Likelihood", y = "Observed maximum heteroplasmy", title = "Neutral variation") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .95),
        legend.key.size = unit(0.5, "cm")) + 
  scale_color_manual(labels = c("G>A and T>C", "All other variants"), values = c("#e41a1c", "#377eb8")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(700, 650), digits = 3, show.legend = FALSE) +
  geom_smooth(method = 'lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))

ggsave("figures/Figure1C.png", width = 5, height = 3)


# SUPPLEMENTARY FIGURE 2B
# Plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD, in tRNA genes only.

ggplot(data = file[grepl("MT-T", file$symbol), ], aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 2) +
  labs(x = "Likelihood", y = "Observed maximum heteroplasmy", title = "Neutral variation in tRNA genes") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .95),
        legend.key.size = unit(0.5, "cm")) + 
  scale_color_manual(labels = c("G>A and T>C", "All other variants"), values = c("#e41a1c", "#377eb8")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(28, 26), digits = 3, show.legend = FALSE) +
  geom_smooth(method='lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))

ggsave("supplementary_figures/FigureS2B.png", width = 5, height = 3)


# SUPPLEMENTARY FIGURE 2C
# Plot the correlation between the mutation likelihood scores and the observed maximum heteroplasmy of neutral variants in gnomAD, in the OriB-OriH region.

ori_file <- read.delim(file = '../output_files/calibration/gnomad_loci_obs_vs_scores_ori.txt', header = TRUE, sep = "\t")

ggplot(data = ori_file, aes(y = obs_max_het, x = sum_likelihood, color = mutation_group)) +
  geom_point(size = 2) +
  labs(x = "Likelihood", y = "Observed maximum heteroplasmy", title = "Neutral variation in OriB-OriH") + 
  paper_theme +
  theme(legend.title = element_blank(),
        legend.position = c(.85, .95),
        legend.key.size = unit(0.5, "cm")) + 
  ylim(c(1, 110)) + 
  scale_color_manual(labels = c("G>A and T>C", "All other variants"), values = c("#e41a1c", "#377eb8")) + 
  stat_cor(method = "pearson", label.x = c(0, 0), label.y = c(110, 100), digits = 3, show.legend = FALSE) +
  geom_smooth(method = 'lm') +
  guides(color = guide_legend(override.aes = list(fill = NA)))

ggsave("supplementary_figures/FigureS2C.png", width = 5, height = 3)


# Save linear model fits, for use in calculating expected.

model_fits <- rbind(c("region", "mutation_group", "item", "value"),
                   c("ref_exc_ori", "G>A_and_T>C", "coefficient",  coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                   c("ref_exc_ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                   c("ref_exc_ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["sum_likelihood"]),
                   c("ref_exc_ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = file[file$mutation_group == "other",]))["(Intercept)"]),
                   c("ori", "G>A_and_T>C", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["sum_likelihood"]),
                   c("ori", "G>A_and_T>C", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "G>A_and_T>C",]))["(Intercept)"]),
                   c("ori", "other", "coefficient", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["sum_likelihood"]),
                   c("ori", "other", "intercept", coef(lm(obs_max_het ~ sum_likelihood, data = ori_file[ori_file$mutation_group == "other",]))["(Intercept)"]))

write.table(model_fits, file = "../output_files/calibration/linear_model_fits.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

