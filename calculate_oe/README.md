## Overview:

`calibrate_model_step1.py`: For each gene/locus, calculate the observed sum maximum heteroplasmy of neutral variants and their sum mutation likelihood scores to fit linear models. User specified argument *-input*, the file with likelihood scores and observed values; the output of `annotate_mutations.py` set as default.

`calibrate_model_step2.R`: Use the output of `calibrate_model_step1.py` to fit linear models for netural variation. Correlation significance testing is also performed to show that the observed level of neutral variation is more significantly correlated with mutation likelihood scores than locus length, using the *cocor* package. 

`oe_functions.py`: Functions used for observed:expected value calculations.

`oe_consequences.py`: Calculate the observed:expected ratios and confidence intervals for different functional classes of variation.

`oe_loci.py`: Calculate the observed:expected ratios and confidence intervals for different functional classes of variation in different genes/loci.

`oe_other.py`: Calculate the observed:expected ratios and confidence intervals for other categories of variation.

`calculate_oe_gnomad.sh`: Apply the model to assess constraint in gnomAD.

`replicate_oe_helix.sh`: Apply the model to a replication dataset, HelixMTdb.

