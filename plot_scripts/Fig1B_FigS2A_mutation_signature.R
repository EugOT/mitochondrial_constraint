library(ggplot2)
library(spgs)

# FIGURE 1B
# plot the mutational signature estimated by the model, across the reference mtDNA (excluding OriB-OriH)

file <- read.delim(file = '../output_files/mutation_likelihoods/mitochondrial_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
# convert to pyridimine mutation
file$pyr_mut <- ifelse(file$mut == "G>T","C>A",
                     ifelse(file$mut == "G>A", "C>T",
                            ifelse(file$mut == "G>C", "C>G",
                                   ifelse(file$mut == "A>T", "T>A",
                                          ifelse(file$mut == "A>C", "T>G",
                                                 ifelse(file$mut == "A>G", "T>C", file$mut))))))
file$pyr_tri <- ifelse(file$REF == "G" | file$REF == "A", reverseComplement(file$trinucleotide, case = "upper"), as.character(file$trinucleotide))
file$strand <- ifelse(file$REF == "G" | file$REF == "A", "Heavy", "Light")
# exclude OriB-OriH m.191-16197
for_plot <- unique(file[file$POS > 191 & file$POS < 16197, c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")]) 
  
ggplot(data = for_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Set2", name = "Strand") + 
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide", y = "Likelihood") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25) 
  
ggsave("figures/Figure1B.png", width = 12, height = 3, dpi = 300)


# SUPPLEMENTARY FIGURE 2A
# plot the mutational signature estimated by the model, across the OriB-OriH region

ori_plot <- unique(file[file$POS < 192 | file$POS > 16196, c("Likelihood", "trinucleotide", "pyr_mut", "pyr_tri", "strand")])

ggplot(data = ori_plot, aes(x = pyr_tri, y = Likelihood)) +
  geom_bar(stat = "identity", aes(fill = strand), position = position_dodge(width = 0.9)) + 
  scale_fill_brewer(palette = "Set2", name = "Strand") + 
  facet_grid(.~pyr_mut, scales = "free") + 
  labs(x = "Trinucleotide", y = "Likelihood in OriB-OriH") + 
  paper_theme + 
  theme(axis.text.x  = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) + 
  geom_vline(xintercept = c(4.525, 8.5, 12.5), linetype = "dashed", colour = "dark grey", size = 0.25) 
  
ggsave("supplementary_figures/FigureS2A.png", width = 12, height = 3, dpi = 300)  

