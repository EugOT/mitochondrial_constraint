library(colorspace)
library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)

# EXTENDED DATA FIGURE 8

# Figure ED8a - stacked bar plot of disease-associated variation by MLC score quartile

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")

# divide into four bins
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))

scores$group <- factor(ifelse(grepl("Cfrm", scores$mitomap_status) | (grepl("athogenic", scores$clinvar_interp)), "Pathogenic",
                              ifelse(grepl("Benign", scores$clinvar_interp), "Benign", "neither")) , 
                       levels = c("Pathogenic", "Benign", "neither"))

# stacked bar plot pathogenic vs benign by score quartile
plotA <- ggplot(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),], 
                aes(group, fill = fct_rev(rank_bin))) + 
  geom_bar(position = "fill", colour = "black") +
  labs(y = "Proportion", fill = "MLC score\nquartile ") + 
  paper_theme +
  theme(axis.title.x = element_blank(),
        plot.margin = unit(c(0.6, 0.15, 0.6, 0.25), "cm")) +
  scale_fill_manual(values = c("#ff4124", "#ffbfaa", "#cfb1ff", "#542eff")) 

# summarize for reference
as.data.frame(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),] 
              %>% group_by(group, rank_bin) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))
as.data.frame(scores[scores$group == "Pathogenic" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence)
                     & !grepl("nr|na", scores$mitomap_plasmy) & scores$mitomap_plasmy != "",] %>% group_by(mitomap_plasmy, rank_bin) %>% summarize(n = n()) %>% mutate(freq = n / sum(n)))


# Figure ED8b - density plot by pathogenic vs benign

plotB <- ggplot(scores[scores$group != "neither" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence),], 
                aes(x = as.numeric(MLC_var_score), fill =  group)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Group') + 
  paper_theme +
  theme(legend.text = element_text(margin = margin(r = 35, unit = "pt")),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "left",
        legend.margin = margin(c(0, 0, 0, 1.5)))


# Figure ED8c - collapse disease plasmy into two groups - plot all pathogenic

scores$plasmy <- ifelse(grepl("\\+\\/", scores$mitomap_plasmy), "At homoplasmy", "Only at heteroplasmy")
plotC <- ggplot(scores[scores$group == "Pathogenic" & grepl("missense|transcript", scores$consequence) & !grepl("gain|lost|terminal", scores$consequence)
                       & !grepl("nr|na", scores$mitomap_plasmy) & scores$mitomap_plasmy != "",], 
                aes(x = as.numeric(MLC_var_score), fill =  plasmy)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Pathogenic with\nMITOMAP plasmy status') + 
  paper_theme +
  xlim(0, 1) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "left",
        legend.margin = margin(c(0, 0, 0, 1.5)))

# write table for curation, for next plots
scores$var <- paste("m.", scores$POS, scores$REF, ">", scores$ALT, sep = "")
#write.table(scores[grepl("Cfrm", scores$mitomap_status), c("var", "symbol", "consequence", "mitomap_status", "mitomap_plasmy", "MLC_var_score")], 
#            file = 'figures/cnfm_by_status.txt', sep = "\t", row.names = FALSE, col.names = TRUE)


# Figure ED8d - collapse disease plasmy into two groups - curated confirmed (above file modified)

curated <- read.delim(file = 'supplementary_datasets/cnfm_by_status_curated.txt', header = TRUE, sep = "\t")
curated$plasmy <- ifelse(grepl("\\+\\/", curated$curated_status), "At homoplasmy", "Only at heteroplasmy")

plotD <- ggplot(curated[grepl("missense|transcript", curated$consequence) & !grepl("gain|lost|terminal", curated$consequence),], aes(x = as.numeric(MLC_var_score), fill = plasmy)) + 
  geom_density(alpha = 0.5) +
  labs(x = 'MLC score', y = 'Density', fill = 'Confirmed pathogenic\nwith curated plasmy status') + 
  paper_theme +
  xlim(0, 1) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.5), "cm"),
        legend.position = "left",
        legend.margin = margin(c(0, 0, 0, 1.5)))


# Figure 8e - boxplot to show the MLC score distribution for indels in population databases

# read in gnomAD, filter to PASS only
gnomad <- read.delim(file = '../required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv', header = TRUE, sep = "\t")
gnomad <- gnomad[gnomad$filters == "PASS",]
gnomad$pos <- gnomad$position
gnomad$type <- ifelse((nchar(as.character(gnomad$ref)) == 1 & nchar(as.character(gnomad$alt)) == 1), "SNV", "indel")

# read in HelixMTdb, extract ref, pos, alt
helix <- read.delim(file = '../required_files/databases/HelixMTdb_20200327.tsv', header = TRUE, sep = "\t")
helix <- helix[str_count(helix$alleles, ",") == 1, ]  # remove multiallelic sites
helix$ref <-  sub(",.*", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$alt <-  sub(".*,", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$pos <- sub("chrM:", "", helix$locus)
helix$type <- ifelse((nchar(as.character(helix$ref)) == 1 & nchar(as.character(helix$alt)) == 1), "SNV", "indel")

# read in MITOMAP, filter to variants observed in genbank
mitomap <- read.delim(file = '../required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt', header = TRUE, sep = "\t")
mitomap <- mitomap[mitomap$gbcnt > 0,]
mitomap$type <- ifelse((nchar(as.character(mitomap$ref)) == 1 & nchar(as.character(mitomap$alt)) == 1) & !grepl(":", mitomap$ref) & !grepl(":", mitomap$alt), "SNV", "indel")

# merge and annotate
gnomad$db <- "gnomAD"
helix$db <- "HelixMTdb"
mitomap$db <- "MITOMAP"
indels <- rbind(gnomad[gnomad$type == "indel", c("pos", "ref", "alt", "db")],
                helix[helix$type == "indel", c("pos", "ref", "alt", "db")],
                mitomap[mitomap$type == "indel", c("pos", "ref", "alt", "db")])

# annotate with position scores
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$pos <- scores$POS
indels <- merge(indels, scores[!duplicated(scores$POS), c("pos", "MLC_pos_score")], by = "pos", all.x = T)

as.data.frame(indels %>% group_by(db) %>% summarize(median = median(MLC_pos_score), n = n()))

plotE <- ggplot(data = indels, aes(x = db, y = MLC_pos_score, fill = db)) + 
  geom_boxplot(position = "dodge") +
  labs(y = 'MLC score', x = "Indels") + 
  paper_theme +
  ylim(0, 1) +
  scale_fill_manual(values = c("#F9F6EE", "#F9F6EE", "#F9F6EE"), guide = FALSE)


# Figure ED8f - plot MLC position scores vs PhyloP

scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
scores$rank_bin <- ifelse(scores$MLC_pos_score <= 0.25, "0.0-0.25", 
                          ifelse(scores$MLC_pos_score > 0.25 & scores$MLC_pos_score <= 0.5, "0.25-0.50", 
                                 ifelse(scores$MLC_pos_score > 0.5 & scores$MLC_pos_score <= 0.75, "0.50-0.75", 
                                        ifelse(scores$MLC_pos_score > 0.75, "0.75-1.0", "error"))))

plotF <- ggplot(data = scores[!duplicated(scores$POS),], aes(x = rank_bin, y = phyloP_score, fill = rank_bin)) + 
  geom_boxplot() +
  labs(x = 'MLC score', y = 'PhyloP score') + 
  paper_theme +
  geom_hline(yintercept = c(0), linetype = "dashed") +
  scale_fill_discrete_sequential(palette = "OrYel", guide = FALSE, order = c(1, 2, 3, 4)) 


# Figure ED8g - the local constraint score distribution of base/amino acids substitutions in gnomAD

# read in and merge in scores - note filter to PASS only in gnomAD
gnomad <- read.delim(file = '../required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv', header = TRUE, sep = "\t")
scores <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
gnomad$var <- paste(gnomad$position, gnomad$ref, gnomad$alt)
scores$var <- paste(scores$POS, scores$REF, scores$ALT)
gnomad <- merge(gnomad[gnomad$filters == "PASS",], scores[, c("var", "MLC_var_score")], by = "var", all.x = T) 

# assign variant type and category
gnomad$type <- ifelse((nchar(as.character(gnomad$ref)) == 1 & nchar(as.character(gnomad$alt)) == 1), "SNV", "indel")
gnomad$category <- factor(ifelse(gnomad$AF_hom == 0, "Heteroplasmy\nonly", 
                                 ifelse(gnomad$AF_hom > 0 & gnomad$AF_hom < 1/50000, "Homoplasmy\nAF <0.002%", 
                                        "Homoplasmy\nAF ≥0.002%")),
                          levels = c("Homoplasmy\nAF ≥0.002%", "Homoplasmy\nAF <0.002%",  "Heteroplasmy\nonly"))
as.data.frame(gnomad %>% group_by(category) %>% summarize(n()))

plotG <- ggplot(data = gnomad[gnomad$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "\ngnomAD SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey", "white"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5))


# Figure ED8h - the local constraint score distribution of base/amino acids substitutions in HelixMTdb

# read in, format and merge in scores
helix <- read.delim(file = '../required_files/databases/HelixMTdb_20200327.tsv', header = TRUE, sep = "\t")
helix <- helix[str_count(helix$alleles, ",") == 1, ]  # remove multiallelic sites
helix$ref <-  sub(",.*", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$alt <-  sub(".*,", "", as.character(gsub('\\[|\\]', '', helix$alleles)))
helix$POS <- sub("chrM:", "", helix$locus)
helix$var <- paste(helix$POS, helix$ref, helix$alt)
helix <- merge(helix, scores[, c("var", "MLC_var_score")], by = "var", all.x = T) 

# assign variant type and category
helix$type <- ifelse((nchar(as.character(helix$ref)) == 1 & nchar(as.character(helix$alt)) == 1), "SNV", "indel")
helix$category <- factor(ifelse(helix$AF_hom == 0, "Heteroplasmy\nonly", 
                                ifelse(helix$AF_hom > 0 & helix$AF_hom < 1/50000, "Homoplasmy\nAF <0.002%", 
                                       "Homoplasmy\nAF ≥0.002%")),
                         levels = c("Homoplasmy\nAF ≥0.002%", "Homoplasmy\nAF <0.002%",  "Heteroplasmy\nonly"))
as.data.frame(helix %>% group_by(category) %>% summarize(n()))

plotH <- ggplot(data = helix[helix$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "\nHelixMTdb SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey", "white"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5))


# Figure ED8i - the local constraint score distribution of base/amino acids substitutions in MITOMAP

# read in, format and merge in scores
mitomap <- read.delim(file = '../required_files/databases/MITOMAP_polymorphisms_2022-07-14.txt', header = TRUE, sep = "\t")
mitomap$var <- paste(mitomap$pos, mitomap$ref, mitomap$alt)
mitomap <- merge(mitomap[!duplicated(mitomap$var),], scores[, c("var", "MLC_var_score")], by = "var", all.x = T)

# assign variant type and category, also calculate allele frequency
mitomap$type <- ifelse((nchar(as.character(mitomap$ref)) == 1 & nchar(as.character(mitomap$alt)) == 1) & !grepl(":", mitomap$ref) & !grepl(":", mitomap$alt), "SNV", "indel")
mitomap$af <- mitomap$gbcnt/56910  # number of sequences in MITOMAP at time of download
mitomap$category <- factor(ifelse(mitomap$af < 1/50000, "AF <0.002%", "AF ≥0.002%"), levels = c("AF ≥0.002%", "AF <0.002%"))
as.data.frame(mitomap %>% group_by(category) %>% summarize(n()))

plotI <- ggplot(data = mitomap[mitomap$type == "SNV",], aes(x = category, y = MLC_var_score, fill = category)) + 
  geom_boxplot(position = "dodge") +
  labs(x = "\nMITOMAP SNVs", y = 'MLC variant score') + 
  paper_theme +
  scale_fill_manual(values = c("dark grey", "light grey"), guide = FALSE) +
  theme(axis.text.x  = element_text(size = 6.5))


# compile panel
ggarrange(
  ggarrange(plotA, 
            ggarrange(plotB, plotC, plotD, nrow = 3, labels = c("b", "c", "d"), font.label = list(size = 10)), 
            ncol = 2, widths = c(0.7, 1), labels = c("a", ""), font.label = list(size = 10)),
  ggarrange(plotE, plotF, ncol = 2, nrow = 1, labels = c("e", "f"), font.label = list(size = 10), widths = c(0.5, 0.5)),
  ggarrange(plotG, plotH, plotI, ncol = 3, nrow = 1, labels = c("g", "h", "i"), widths = c(1, 1, 0.75), font.label = list(size = 10)),
  nrow = 3, ncol = 1, heights = c(0.8, 0.45, 0.55))

ggsave("extended_data_figures/FigureED8.jpeg", width = 180, height = 170, dpi = 600, units = c("mm"))


