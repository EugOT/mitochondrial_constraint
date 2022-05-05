## Overview:

`compile_denovo.py`: Parse lists of de novo variants obtained from the literature, and in-house datasets. Details on each data source are provided in <placeholder link to manuscript>. Outputs list of all de novo and their source.

`compare_denovo.py`: Compare the mutational likelihoods of transitions, across different sources (germline, somatic tissue, somatic cancer) and sample de novo counts.

`filter_denovo.py`: Remove any de novo from outlier samples. Outputs list of de novo mutations used to calculate mutability. User specified arguments *-germline_max*, *-som_tissue_max*, *-som_cancer_max* to indicate the maximum sample de novo count used for filtering for each source. Defaults provided based on analyses.

`composite_likelihood_mito.py`: Code to apply the mitochondrial composite likelihood model. Outputs file with mutational likelihood scores for every single nucleotide variant in the mtDNA. User specified arguments *-context_size* to indicate how many nucleotides to include for sequence context, 3 for trinucleotide set as default, and *-denovo_list* the path to the list of de novo to use, the output of `filter_denovo.py` set as default.

`annotate_mutations.py`: Annotate the output of `composite_likelihood_mito.py` with annotations needed for downstream analyses, such as variant consequence, gene, and in silico predictions. User specified argument *-input*; the output of `composite_likelihood_mito.py` set as default.

`calibrate_model.py`: Fit linear models of the observed sum maximum heteroplasmy of neutral variants and their mutation likelihood scores. User specified argument *-input*, the file with likelihood scores and observed maximum heteroplasmy; the output of `annotate_mutations.py` set as default.

`correlation_significance_testing.R`: Show that the observed level of neutral variation is more significantly correlated with mutation likelihood scores than locus length, using the *cocor* package. 

`simulate_heteroplasmy.R`: Apply a computational model of germline mtDNA mutation and heteroplasmy drift to support a correlation between mutation rates and maximum heteroplasmy.
