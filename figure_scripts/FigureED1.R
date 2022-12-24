library(forcats)
library(ggplot2)
library(ggpubr)
library(spgs)
library(stringr)
library(tidyverse)

# EXTENDED DATA FIGURE 1

# Figure ED1a - plot the mutational signature estimated by the model, across the OriB-OriH region

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# convert to pyridimine mutation
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
file$pyr_mut <- ifelse(file$mut == "G>T","C>A",
                     ifelse(file$mut == "G>A", "C>T",
                            ifelse(file$mut == "G>C", "C>G",
                                   ifelse(file$mut == "A>T", "T>A",
                                          ifelse(file$mut == "A>C", "T>G",
                                                 ifelse(file$mut == "A>G", "T>C", file$mut))))))
file$pyr_tri <- ifelse(file$REF == "G" | file$REF == "A", reverseComplement(file$trinucleotide, case = "upper"), as.character(file$trinucleotide))
file$strand <- factor(ifelse(file$REF == "G" | file$REF == "A", "Heavy", "Light"), levels = c("Light", "Heavy"), labels = c("Reference / Light", "Reverse complement / Heavy"))

# subset to ori
ori_plot <- unique(file[file$POS < 192 | file$POS > 16196, c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")])

plotA <- ggplot(data = ori_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Pastel1", name = "Strand") +
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide in OriB-OriH", y = "Likelihood") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(3, "mm")) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25)


# Figure ED1b-c - disease-associated variants by consequence 

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$consequence_two <- factor(ifelse(grepl("MT-T", file$symbol), 'tRNA', 
                                      ifelse(grepl("MT-R", file$symbol), 'rRNA', 
                                             ifelse(grepl(",|lost|retained|incomplete|intergenic", file$consequence), "Other", 
                                                    as.character(file$consequence)))),
                               levels = c("synonymous_variant", "missense_variant", "stop_gained", "tRNA", "rRNA", "Other"),
                               labels = c("Synonymous", "Missense", "Stop gain", "tRNA", "rRNA", "Other"))
clinvar_included <- c("Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic")
plot_data_clinvar <- as.data.frame(file[file$clinvar_interp %in% clinvar_included,] %>% group_by(consequence_two) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
plot_data_mitomap <- as.data.frame(file[grepl("Cfrm|Reported", file$mitomap_status),] %>% group_by(consequence_two) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))

mycolors = c('#4daf4a', '#377eb8', '#b22222', '#ff7f00', '#984ea3', 'grey')

plotB <- ggplot(plot_data_clinvar, aes(x = "", y = n, fill = consequence_two)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycolors) +
  labs(fill = "ClinVar\nconsequences") +
  geom_text(aes(x = 1.6, label = ifelse(consequence_two == "rRNA", "      3%", ifelse(consequence_two == "Synonymous", "5%   ", paste0(round(freq * 100, 0), "%")))), 
            position = position_stack(vjust = .5), size = 2.5) +
  theme_void() +
  theme(plot.margin = unit(c(0.25, 0.15, 0.25, 0.15), "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(4.5, "mm"))

plotC <- ggplot(plot_data_mitomap, aes(x = "", y = n, fill = consequence_two)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = mycolors) +
  labs(fill = "MITOMAP\nconsequences") +
  geom_text(aes(x = 1.6, label = paste0(round(freq * 100, 0), "%")), 
            position = position_stack(vjust = .5), size = 2.5) +
  theme_void() +
  theme(plot.margin = unit(c(0.25, 0.15, 0.25, 0.15), "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(4.5, "mm"))


# Figure ED1d - plot the observed:expected ratio and 90% confidence interval for different categories of in silico prediction

file <- read.delim(file = '../output_files/oe/insilicos_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file$group <- factor(str_split(file$prediction, "\\-", simplify = T)[,2], 
                     levels = c("APOGEE", "MitoTip", "HmtVar")) 
file$group2 <- ifelse(file$group == "APOGEE", "Missense", "tRNA")
file$prediction <- factor(str_split(file$prediction, "\\-", simplify = T)[,1], 
                          levels=c("polymorphic", "likely_polymorphic", "likely_pathogenic", "pathogenic", "likely benign", "possibly benign", "possibly pathogenic", "likely pathogenic", "Neutral", "Pathogenic"),
                          labels=c("Polymorphic", "Likely polymorphic", "Likely pathogenic", "Pathogenic", "Likely benign", "Possibly benign", "Possibly pathogenic", "Likely pathogenic", "Neutral", " Pathogenic"))

plotD <- ggplot(file, aes(y = prediction, x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.9, colour = "#777575") +
  geom_point(aes(colour = group2), shape = 18, size = 6) +
  labs(x = "observed:expected ratio") + 
  paper_theme +
  facet_grid(rows = vars(group), scales = "free", space = 'free') +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.y = element_text(size = 7)) +
  scale_color_manual(values = c('#377eb8', '#ff7f00'), name = "Algorithm for:") +
  guides(color = FALSE) 


# Figure ED1e - plot the observed:expected ratio and 90% confidence interval for functional classes of variation in mtDNA in HelixMTdb

file <- read.delim(file = '../output_files/oe/replication_dataset/helix_consequences_obs_exp.txt', header = TRUE, sep = "\t", stringsAsFactors = FALSE)
file <- file[!grepl("lost", file$consequence), ]

file$group <- factor(ifelse(grepl("RNA", file$consequence), "RNA", 
                            ifelse(file$consequence == "intergenic", "Noncoding", 
                                   "Protein")),
                     levels = c("Protein", "RNA", "Noncoding"),
                     labels = c("Protein", "RNA", "Other"))
file$consequence <- factor(file$consequence, levels = c("synonymous", "missense", "stop_gain", "tRNA", "rRNA", "intergenic"), 
                           labels = c("Synonymous", "Missense", "Stop gain", "tRNA", "rRNA", "Intergenic"))

mycolors = c('#4daf4a', '#377eb8', '#b22222', '#ff7f00', '#984ea3', '#ffcc00')

plotE <- ggplot(file, aes(y = fct_rev(consequence), x = as.numeric(obs.exp))) + 
  geom_errorbar(aes(xmin = as.numeric(lower_CI), xmax = as.numeric(upper_CI)), width = 0.65, colour = "#777575") +
  geom_point(aes(colour = consequence), shape = 18, size = 6, show.legend = FALSE) +
  labs(x = "observed:expected ratio in HelixMTdb") + 
  paper_theme +
  theme(axis.title.y = element_blank(),
        panel.background = element_rect(fill = NA, color = "grey", linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  facet_grid(rows = vars(group), scales="free", space='free') + 
  scale_color_manual(values = mycolors)


# compile panel
ggarrange(plotA, 
          ggarrange(plotB, plotC, nrow = 1, ncol = 2, widths = c(1, 1), labels = c("b", "c"), font.label = list(size = 10)),
          ggarrange(plotD, plotE, nrow = 1, ncol = 2, widths = c(1, 0.9), labels = c("d", "e"), font.label = list(size = 10)),
          nrow = 3, ncol = 1, heights = c(0.75, 0.65 , 1), labels = c("a", "", ""), font.label = list(size = 10))

ggsave("extended_data_figures/FigureED1.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))  



