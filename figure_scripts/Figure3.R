library(forcats)
library(ggplot2)
library(ggpubr)
library(png)
library(stringr)

# FIGURE 3

# Figure 3a - plot the observed:expected ratio and 90% confidence interval for variation in each RNA

file <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep = "\t")
file <- file[grepl("MT-T|MT-R", file$locus),]
file$consequence <- factor(ifelse(grepl("tRNA", file$description), "tRNA", 
                                  ifelse(grepl("MT-R", file$locus), "rRNA", 
                                         as.character(file$consequence))), levels = c("tRNA", "rRNA"))

# use RNA colors
mycolors = c('#ff7f00', '#984ea3')

plotA <- ggplot(file, aes(y = reorder(locus, -upper_CI), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.95, colour = "#777575") +
  geom_point(aes(colour = consequence), shape = 18, size = 4.5) +
  labs(x = "Observed:expected ratio") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_text(size = 7),
        plot.margin = unit(c(0.25, 0.25, 0.15, 0.25), "cm")) +
  xlim(0, 1) +
  scale_color_manual(values = mycolors, name = "RNA type") +
  facet_grid(rows = vars(consequence), scales = "free", space = 'free') +
  guides(colour = FALSE)


# Figure 3b - plot the observed:expected ratio and 90% CI for RNA base types

file <- read.delim(file = '../output_files/oe/RNA_base_types_obs_exp.txt', header = TRUE, sep = "\t")
file$group <- factor(ifelse(grepl("tRNA", file$RNA_base_type), "tRNA", "rRNA"), levels = c("tRNA", "rRNA"))
file$RNA_base_type <- factor(str_sub(file$RNA_base_type, start = 1, end = -6), 
                             levels = c("WC", "non-WC", "loop-or-other", "WC", "non-WC", "loop-or-other", "modified", "modified", "non-modified", "non-modified"), 
                             labels = c("WC pair", "Non-WC pair", "Loop or other", "WC pair", "Non-WC pair", "Loop or other", "Modified", "Modified", "Non-modified", "Non-modified"))
# use RNA colors
mycolors = c('#ff7f00', '#984ea3')

# base types
plotB <- ggplot(file[!grepl("odified", file$RNA_base_type),], aes(y = fct_rev(RNA_base_type), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.9, colour = "#777575") +
  geom_point(aes(colour = group), shape = 18, size = 5.5) +
  labs(x = "Observed:expected ratio") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.15, 0.25), "cm")) +
  xlim(0, 1) +
  scale_color_manual(values = mycolors, name = "RNA type") +
  facet_grid(rows = vars(group), scales = "free", space = 'free') +
  guides(colour = FALSE)


# Figure 3c - plot the observed:expected ratio and 90% CI for tRNA domains

file <- read.delim(file='../output_files/oe/tRNA_domains_obs_exp.txt', header = TRUE, sep = "\t")
file$tRNA_domain <- factor(file$tRNA_domain, levels = c("", "Acceptorstem", "D-stem", "D-loop", "Anticodonstem", "Anticodonloop", "Variableregion", "T-stem", "T-loop"), 
                           labels = c("", "Acceptor stem", "D-stem", "D-loop", "Anticodon stem", "Anticodon loop", "Variable region", "T-stem", "T-loop"))
file$group <- "tRNA"

plotC <- ggplot(file[file$tRNA_domain != "",], aes(y = fct_rev(tRNA_domain), x = as.numeric(obs.exp), color = tRNA_domain)) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.85, colour = "#777575") +
  geom_point(shape = 18, size = 5.5) +
  labs(x = "Observed:expected ratio") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        plot.margin = unit(c(0.35, 0.25, 0.15, 0.0), "cm")) +
  xlim(0, 1) +
  facet_grid(rows = vars(group), scales = "free", space = 'free') +
  guides(colour = FALSE)


# Figure 3d-e - tRNA positional constraint in 2D and 3D structures

svg <- readPNG("figures/Figure3d_tRNA_pos.png")  # 2D
plotD <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.25, 0.0, 0.15, 0.1), "cm"))

svg <- readPNG("figures/Figure3e_tRNA_3D.png")  # 3D
plotE <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.75, 0.0, 0.25, 0.2), "cm"))


# Figure 3f - regional constraint in MT-RNR1

svg <- readPNG("figures/Figure3f_RNR1.png") 
plotF <- ggplot() + 
  background_image(svg) +
  theme(plot.margin = unit(c(0.15, 0.05, 0.15, 0.0), "cm"))


# assemble panel 
ggarrange(
  ggarrange(plotA,
            ggarrange(plotB, plotC, labels = c("b", "c"), nrow = 2, ncol = 1, heights = c(0.45, 0.55), font.label = list(size = 10)),
            nrow = 1, ncol = 2, labels = c("a", ""), font.label = list(size = 10)),
  ggarrange(plotD, NULL, plotE, NULL, plotF, labels = c("d", "", "e", "", "f"), nrow = 1, ncol = 5, font.label = list(size = 10), widths = c(0.55, -0.01, 0.35, -0.015, 0.5)),
  nrow = 2, ncol = 1, heights = c(1.05, 1))

ggsave("figures/Figure3.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))

