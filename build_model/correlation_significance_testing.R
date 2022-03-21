library(cocor)

# Show that the observed level of neutral variation is more significantly correlated with mutation likelihood scores than locus length, using the cocor package. 
# Testing two correlations which have dependent groups (ie they have one variable in common), and a one-sided hypothesis
# cor.test is used to extract the correlation coefficient for testing

file <- read.delim(file = '../output_files/calibration/gnomad_loci_obs_vs_scores.txt', header = TRUE, sep = "\t")

sink("../output_files/calibration/correlation_significance_testing.txt", append = FALSE)
print("G>A and T>C mutation group (for reference excluding OriB-OriH)")
cocor.dep.groups.overlap(cor.test(file[file$mutation_group == "G>A_and_T>C", c("obs_max_het")], file[file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], method = "pearson")$estimate, 
                         cor.test(file[file$mutation_group == "G>A_and_T>C", c("obs_max_het")], file[file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                         cor.test(file[file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], file[file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                         n = 39, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0)

print("All other mutation group (for reference excluding OriB-OriH)")
cocor.dep.groups.overlap(cor.test(file[file$mutation_group == "other", c("obs_max_het")], file[file$mutation_group == "other", c("sum_likelihood")], method = "pearson")$estimate, 
                         cor.test(file[file$mutation_group == "other", c("obs_max_het")], file[file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                         cor.test(file[file$mutation_group == "other", c("sum_likelihood")], file[file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                         n = 39, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0)

ori_file <- read.delim(file = '../output_files/calibration/gnomad_loci_obs_vs_scores_ori.txt', header = TRUE, sep = "\t")

print("G>A and T>C mutation group for OriB-OriH region")
cocor.dep.groups.overlap(cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("obs_max_het")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], method = "pearson")$estimate, 
                         cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("obs_max_het")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                         cor.test(ori_file[ori_file$mutation_group == "G>A_and_T>C", c("sum_likelihood")], ori_file[ori_file$mutation_group == "G>A_and_T>C", c("length")], method = "pearson")$estimate, 
                         n = 8, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0)

print("All other mutation group for OriB-OriH region")
cocor.dep.groups.overlap(cor.test(ori_file[ori_file$mutation_group == "other", c("obs_max_het")], ori_file[ori_file$mutation_group == "other", c("sum_likelihood")], method = "pearson")$estimate, 
                         cor.test(ori_file[ori_file$mutation_group == "other", c("obs_max_het")], ori_file[ori_file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                         cor.test(ori_file[ori_file$mutation_group == "other", c("sum_likelihood")], ori_file[ori_file$mutation_group == "other", c("length")], method = "pearson")$estimate, 
                         n = 8, alternative = "greater", alpha = 0.05, conf.level = 0.95, null.value = 0)
sink()