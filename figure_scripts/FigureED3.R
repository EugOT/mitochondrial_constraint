library(ggplot2)
library(ggpubr)
library(ggrepel)

# EXTENDED DATA FIGURE 3

# Figure ED3a-b - plot the mean phyloP score against OEUF for each protein (a) and RNA (b) gene

# mean phyloP scores
file <- read.delim(file = '../output_files/other/gene_phylop.txt', header = TRUE, sep = "\t")
# merge with OEUF values
gene_oe <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep = "\t")
gene_oe <- merge(gene_oe[gene_oe$consequence == "missense" | gene_oe$consequence == "SNV", c("locus", "upper_CI")], file, by = "locus", all.x = T) 

# label theme for these plots
paper_repel <- geom_text_repel(aes(label = locus), 
                               color = 'grey50',
                               size = 1.75,
                               box.padding = 0.5,
                               point.padding = 0.5,
                               segment.color = 'grey',
                               segment.size = 0.25, 
                               nudge_x = -0.05)

# protein genes
plotA <- ggplot(gene_oe[grepl("MT-A|MT-C|MT-N", gene_oe$locus),], aes(x = upper_CI, y = median_phylop_score)) + 
  geom_point(color = '#377eb8') +
  labs(y = "Median phyloP score", x = "Missense OEUF") + 
  paper_theme +
  xlim(0, 1) +
  ylim(-2, 8) +
  paper_repel +
  stat_cor(method = "pearson", digits = 2, show.legend = FALSE, size = 2.5)

# RNA genes
gene_oe$type <- ifelse(grepl("MT-T", gene_oe$locus), "tRNA", "rRNA")
plotB <- ggplot(gene_oe[grepl("MT-T|MT-R", gene_oe$locus),], aes(x = upper_CI, y = median_phylop_score)) + 
  stat_cor(aes(label = locus), method = "pearson", digits = 2, show.legend = FALSE, size = 2.5) +
  geom_point(aes(color = type)) +
  labs(y = "Median phyloP score", x = "OEUF") + 
  xlim(0, 1) +
  ylim(-2, 8) +
  paper_theme +
  paper_repel +
  scale_color_manual(values = c('#984ea3', '#ff7f00')) +
  theme(legend.position = "none") 
  

# Figure ED3c - plot each tRNA's codon usage against their OEUF

# codon usage
file <- read.delim(file = '../output_files/other/tRNA_codon_counts.txt', header = TRUE, sep = "\t")
file$locus <- file$tRNA_codon
# merge with OEUF values
gene_oe <- read.delim(file='../output_files/oe/genes_obs_exp.txt', header = TRUE, sep = "\t")
gene_oe <- merge(gene_oe[grepl("MT-T", gene_oe$locus), c("locus", "upper_CI")], file, by = "locus", all.x = T) 

plotC <- ggplot(gene_oe, aes(x = upper_CI, y = count)) + 
  geom_point(color = '#377eb8') +
  labs(y = "tRNA codon count", x = "OEUF") + 
  xlim(0, 1) +
  paper_theme +
  paper_repel +
  stat_cor(method = "pearson", digits = 2, show.legend = FALSE, size = 2.5) +
  theme(plot.margin = unit(c(0.25, 3, 0.25, 3), "cm"))


# compile panel
ggarrange(
  ggarrange(plotA, plotB, labels = c("a", "b"), font.label = list(size = 10)),
  plotC, nrow = 2, ncol = 1, labels = c("", "c"), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED3.jpeg", width = 180, height = 150, dpi = 600, units = c("mm"))



  