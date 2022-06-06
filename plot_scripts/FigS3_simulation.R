library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)


# SUPPLEMENTARY FIGURE 3

# Figure S3A - violin plot to show distribution of heteroplasmy levels after one generation, to replicate Figure 2B for mature oocytes from Colgnahi et al

results <- read.table("../output_files/simulation/colnaghi_etal_replication.txt", header = TRUE)

plotA <- ggplot(data = results, aes(y = heteroplasmy, x = fct_rev(as.factor(bottleneck_size)), fill = as.factor(bottleneck_size))) +
  geom_violin(scale = "width") + 
  labs(x = "Bottleneck size", y = "Heteroplasmy level") +
  paper_theme +
  guides(fill = FALSE) +
  stat_summary(fun = "mean", geom = "point", color = "red") +
  scale_fill_brewer(palette = "Blues")


# Figure S3B - plotting the change in heteroplasmy between generations, as a replication of Figure 2D from Wei et al
# will show for midpoint mutation rate, 1e-08
# note Wei et al define de novo as being >1% heteorplasmy in offspring and not detected in mother, will apply same criteria

results <- read.table("../output_files/simulation/simulation_results.txt", header = TRUE)

# create a dataframe with heteroplasmy values for the same lineage across two subsequent generations
table1 <- results[results$generation == 4 & results$mutation_rate == 1e-08,]  # ie mother
table2 <- results[results$generation == 5 & results$mutation_rate == 1e-08,]  # ie offspring
colnames(table2) <- c("individual2", "heteroplasmy2", "generation2", "mutation_rate2", "bottleneck_size2", "back_mutation_rate2", "starting_heteroplasmy2")
joint_table <- cbind(table1[order(table1$individual),], table2[order(table2$individual),])
joint_table$change <- joint_table$heteroplasmy2 - joint_table$heteroplasmy

# remove any homoplasmic variants that are fixed
joint_table <- joint_table[joint_table$heteroplasmy != 1 & joint_table$heteroplasmy2 != 1,]

# label for plotting
joint_table$label <- ifelse(joint_table$heteroplasmy < 0.01 & joint_table$heteroplasmy2 > 0.01, "de novo",
                          ifelse(joint_table$heteroplasmy > 0.01 & joint_table$heteroplasmy2 < 0.01,"lost",
                                 ifelse(joint_table$heteroplasmy > 0.01 & joint_table$heteroplasmy2 > 0.01,"transmitted",
                                        ifelse(joint_table$heteroplasmy < 0.01 & joint_table$heteroplasmy2 < 0.01, "no mut", "error"))))

plotB <- ggplot(data = joint_table[joint_table$label != "no mut", ], aes(color = label, x = reorder(individual, change), change)) +
  geom_point(stat = "identity", aes(shape = 124)) +
  labs(x = "Heteroplasmic variant rank", y = "Heteroplasmy change") + 
  paper_theme + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous(limits = c(-0.75, 1), breaks = seq(-0.75, 1, by = 0.25)) +
  scale_shape_identity()


# Figure S3C - plotting the maximum heteroplasmy in all 10,000 lineages across generations

max_het <- aggregate(results$heteroplasmy, by = list(results$generation, results$mutation_rate), max)
max_het <- rbind(max_het, 
                 c(0, "1e-07", 0), c(0, "7.5e-08", 0), c(0, "5e-08", 0), c(0, "2.5e-08", 0), c(0, "1e-08", 0), c(0, "7.5e-09", 0), c(0, "5e-09", 0), c(0, "2.5e-09", 0), c(0, "1e-09", 0))  # add 0 for generation 0
max_het$Group.2 <- factor(max_het$Group.2, levels = c("1e-07", "7.5e-08", "5e-08", "2.5e-08", "1e-08", "7.5e-09", "5e-09", "2.5e-09", "1e-09"))

plotC <- ggplot(data = max_het, aes(y = as.numeric(x), x = as.numeric(Group.1), color = fct_rev(as.factor(Group.2)))) +
  geom_point() +
  geom_line() + 
  labs(x = "Generation", y = "Maximum heteroplasmy in population") + 
  paper_theme +
  ylim(c(0, 1.0)) +
  labs(color = "Mutation rate")


# Figure S3D - plot bootstrap data to show correlation between mutation rate and maximum heteroplasmy

results <- read.table("../output_files/simulation/sampling_max_heteroplasmy.txt", header = TRUE)
results$mutation_rate <- factor(results$mutation_rate, levels = c("1e-07", "7.5e-08", "5e-08", "2.5e-08", "1e-08", "7.5e-09", "5e-09", "2.5e-09", "1e-09"))

plotD <- ggplot(data = results, aes(y = as.numeric(max_heteroplasmy), x = fct_rev(mutation_rate), fill = fct_rev(mutation_rate))) +
  geom_boxplot() + 
  labs(x = "Mutation rate", y = "Maximum heteroplasmy") + 
  paper_theme +
  guides(fill = FALSE)


# collate figure 
ggarrange(plotA, plotB, plotC, plotD, ncol = 2, nrow = 2, heights = c(1, 1.15), labels = c('a', 'b', 'c', 'd'))
ggsave("supplementary_figures/FigureS3.png", width = 10, height = 6, dpi = 300)

