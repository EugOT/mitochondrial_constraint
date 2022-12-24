# Supplementary Dataset 1 - mutational likelihoods

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
file$mut <- paste(file$REF, ">", file$ALT, sep = "")
file$Likelihood <- round(file$Likelihood, 3)

# write references excluding ori and ori region separately
file$group <- "reference_excluding_ori"
write.table(unique(file[file$POS > 191 & file$POS < 16197, c("trinucleotide", "mut", "Likelihood", "group")]), 
            file = 'supplementary_datasets/supplementary_dataset_1.txt', col.names = c("Trinucleotide", "Mutation", "Likelihood", "Reference"), row.names = FALSE, sep = '\t', append = FALSE)
file$group <- "ori_region"
write.table(unique(file[file$POS < 192 | file$POS > 16196, c("trinucleotide", "mut", "Likelihood", "group")]), 
            file = 'supplementary_datasets/supplementary_dataset_1.txt', col.names = FALSE, row.names = FALSE, sep = '\t', append = TRUE)


# Supplementary Dataset 2 - de novo mutations used to build model

file <- read.delim(file = '../output_files/denovo/all_denovo.txt', header = TRUE, sep= "\t")

# will provide the list that includes sample details, this was filtered (1) - only include cancer samples with 1 de novo
file$exclude <- ifelse(grepl("cancer", file$sample) & file$sample_denovo_count > 1, "yes", "no")

# and (2) to only include base substitutions in reference sequence (i.e. remove any annotation errors)
snvs <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")
snvs$var <- paste(snvs$REF, snvs$POS, snvs$ALT, sep = "")
file$snv <- ifelse(file$denovo %in% snvs$var, "yes", "no")

# write filtered list to file
write.table(file[file$exclude == "no" & file$snv == "yes", c("denovo", "sample")], 
            file = 'supplementary_datasets/supplementary_dataset_2.txt', col.names = c("Denovo", "Sample"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 3 - constraint metrics for each gene

file <- read.delim(file = '../output_files/oe/genes_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
file$consequence <- ifelse(file$consequence == "SNV", "RNA_variant", as.character(file$consequence))
write.table(file[, c("locus", "start", "end", "consequence", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_3.txt', col.names = c("Symbol", "Start_position", "End_position", "Consequence", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 4 - regional constraint intervals

file <- read.delim(file = '../output_files/regional_constraint/final_regional_constraint_intervals.txt', header = TRUE, sep= "\t")
file$obs_max_het <- round(file$obs_max_het, 3)
file$exp_max_het <- round(file$exp_max_het, 3)
file$ratio_oe <- round(file$ratio_oe, 3)
write.table(file[, c("locus", "start", "end", "protein_pos_start", "protein_pos_end", "obs_max_het", "exp_max_het", "ratio_oe", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_4.txt', col.names = c("Symbol", "Start_position", "End_position", "Protein_residue_start", "Protein_residue_end", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 5 - curated missense from VCGS

mcri <- read.delim(file='../../../patient_VUS/MCRI/mtDNA_variants_combined_final_copy.txt', header = TRUE, sep = "\t")
mcri$var <- paste(mcri$REF, mcri$POS, mcri$ALT)
mcri$Classification <- factor(ifelse(grepl("Class 5|Class 4", mcri$Classification), "Pathogenic & Likely pathogenic",
                            ifelse(grepl("Class 1|Class 2", mcri$Classification), "Benign & Likely Benign", paste("VUS ", as.character(mcri$Classification), sep = ""))),
                     levels = c("Benign & Likely Benign", "VUS Class 3c", "VUS Class 3b", "VUS Class 3a", "Pathogenic & Likely pathogenic"),
                     labels = c("Benign & Likely Benign", "VUS of low clinical significance", "VUS ", 
                                "VUS of high clinical significance", "Pathogenic & Likely pathogenic"))
# annotate with consequences
file <- read.delim(file='../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
file$var <- paste(file$REF, file$POS, file$ALT)

# restrict to missense, as most severe consequence 
mcri <- merge(mcri, file[, c("var", "consequence")], by = "var", all.x = T) 
mcri <- mcri[grepl("missense", mcri$consequence) & !grepl("stop_gain|start_lost|stop_lost|incomplete_terminal", mcri$consequence),]

write.table(mcri[order(mcri$POS), c("POS", "REF", "ALT", "Classification")], 
            file = 'supplementary_datasets/supplementary_dataset_5.txt', col.names = c("Position", "Reference", "Alternate", "Classification_group"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 6 - Constraint metrics for each position within the tRNA secondary structure

file <- read.delim(file = '../output_files/oe/tRNA_position_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
# manually reorder
write.table(file[c(2:23,75,24:48,76,49:52,71:72,53:59,73:74,60:70), c("tRNA_position", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_6.txt', col.names = c("tRNA_position", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 7 - Constraint metrics for non-coding elements

file <- read.delim(file = '../output_files/oe/noncoding_obs_exp.txt', header = TRUE, sep= "\t")
file$observed <- round(file$observed, 3)
file$expected <- round(file$expected, 3)
file$obs.exp <- round(file$obs.exp, 3)
write.table(file[, c("locus", "description", "start", "end", "observed", "expected", "obs.exp", "lower_CI", "upper_CI")], 
            file = 'supplementary_datasets/supplementary_dataset_7.txt', col.names = c("Locus", "Description", "Start_position", "End_position", "Observed", "Expected", "obs:exp", "Lower_CI", "Upper_CI"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 8 - MLC scores for each base position

file <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
file$MLC_pos_score <- round(file$MLC_pos_score, 4)
write.table(unique(file[, c("POS", "MLC_pos_score")]), 
            file = 'supplementary_datasets/supplementary_dataset_8.txt', col.names = c("Position", "MLC_pos_score"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 9 - MLC scores for every SNV

file <- read.delim(file = '../output_files/local_constraint/per_base_local_constraint.txt', header = TRUE, sep = "\t")
file$MLC_var_score <- round(file$MLC_var_score, 4)
write.table(file[, c("POS", "REF", "ALT", "consequence", "MLC_var_score")], 
            file = 'supplementary_datasets/supplementary_dataset_9.txt', col.names = c("Position", "Reference", "Alternate", "Consequence", "MLC_score"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 10 - curated confirmed MITOMAP variants by plasmy

file <- read.delim(file = 'supplementary_datasets/cnfm_by_status_curated.txt', header = TRUE, sep = "\t")
file$mitomap_plasmy <- paste("''", file$mitomap_plasmy, "'", sep = "")  # need '' for google sheets
file$curated_status <- paste("''", file$curated_status, "'", sep = "")
file$curated_homoplasmy_report <- gsub(" ", "", file$curated_homoplasmy_report)
write.table(file[grepl("missense|transcript", file$consequence) & !grepl("gain|lost|terminal", file$consequence), c("var", "symbol", "mitomap_plasmy", "curated_status", "curated_homoplasmy_report")], 
            file = 'supplementary_datasets/supplementary_dataset_10.txt', col.names = c("Variant", "Symbol", "MITOMAP_plasmy", "Curated_status", "Curated_homoplasmy_report"), row.names = FALSE, sep = '\t')


# Supplementary Dataset 11 - all possible SNVs and their annotations used in manuscript

file <- read.delim(file = '../output_files/mutation_likelihoods/mito_mutation_likelihoods_annotated.txt', header = TRUE, sep = "\t")

# merge in regional constraint annotations
file$var <- paste(file$REF, file$POS, file$ALT)
rc <- read.delim(file='../output_files/regional_constraint/mito_regional_constraint_annotation.txt', header = TRUE, sep = "\t")
rc$var <- paste(rc$REF, rc$POS, rc$ALT)
file <- merge(file, rc[, c("var", "in_rc", "min_distance_to_rc")], by = "var", all.x = T) 
file <- file[order(file$POS),]

write.table(file[, c("POS", "REF", "ALT", "symbol", "consequence", "protein_position", "gnomad_max_hl", "gnomad_af_hom", 
                     "in_phylotree", "phyloP_score", "tRNA_position", "tRNA_domain", "RNA_base_type", "RNA_modified", "rRNA_bridge_base", "uniprot_annotation", "other_prot_annotation", 
                     "apogee_class", "mitotip_class", "hmtvar_class", 
                     "helix_max_hl", "helix_af_hom", "mitomap_af", 
                     "mitomap_status", "mitomap_plasmy", "clinvar_interp", 
                     "in_rc", "min_distance_to_rc")], 
            file = 'supplementary_datasets/supplementary_dataset_11.txt', 
            col.names = c("Position", "Reference", "Alternate", "Symbol", "Consequence", "Protein_position", 
                          "gnomad_max_hl", "gnomad_af_hom", 
                          "in_phylotree", "phyloP_score", "tRNA_position", "tRNA_domain", "RNA_base_type", "RNA_modified", "rRNA_bridge_base", "uniprot_annotation", "other_protein_annotation", 
                          "apogee_class", "mitotip_class", "hmtvar_class", 
                          "helix_max_hl", "helix_af_hom", "mitomap_af", 
                          "mitomap_status", "mitomap_plasmy", "clinvar_interp", 
                          "in_rc", "min_distance_to_rc"), row.names = FALSE, sep = '\t')

