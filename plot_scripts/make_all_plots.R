# generate all figures
dir.create("figures")
dir.create("supplementary_figures")

library(ggplot2)

# ggplot theme to use for manuscript
paper_theme <- theme(axis.title.x = element_text(size = 10),
                     axis.text.x  = element_text(size = 10),
                     axis.title.y = element_text(size = 10),
                     axis.text.y  = element_text(size = 10),
                     axis.ticks.length = unit(0.15, "cm"),
                     panel.background = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     legend.title = element_text(size = 10),
                     legend.text = element_text(size = 10))

# main text figures
source("Fig1B_FigS2A_mutation_signature.R", echo = TRUE) 
source("Fig1C_FigS2B-C_calibrate_model.R", echo = TRUE) 

# supplementary figures
source("FigS1_compare_denovo.R", echo = TRUE) 
