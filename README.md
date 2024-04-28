# Code for mitochondrial genome constraint manuscript

This repository provides the scripts used to quantify constraint across the human mitochondrial genome, and to generate the data and figures presented in [our manuscript.](<https://www.biorxiv.org/content/10.1101/2022.12.16.520778v2>) We assessed constraint by identifying genes and regions where the observed mtDNA variation in gnomAD is less than expected, under neutrality.

Note that pre-computed mitochondrial constraint metrics are available in the Supplementary Dataset of the manuscript. 

## Overview:

The workflow used is reflected below. See the README in each directory for more information on contents.

`required_files`: input files for code.

`build_model`: code used to build the mitochondrial mutational model, and to curate of de novo mutations used to quantify mutability.

`calculate_oe`: code for assessment of mitochondrial constraint across functional variant classes and loci. 

`regional_constraint`: code used to identify regional constraint.

`local_constraint`: code used to identify local constraint and generate the mitochondrial local constraint (MLC) score.

`other`: code for other analyses presented in the manuscript.

`figure_scripts`: scripts used for generating figures, video and datasets.

The code was run using Python v3.10, R v3.6.1, and ChimeraX v1.3 on a MacOS and or Linux system.

A copy of this repository can be obtained on a local computer by downloading as a zip file or cloning using Git.

## Reference

[1] Nicole J. Lake, Wei Liu, Stephanie L. Battle, Kristen M. Laricchia, Grace Tiao, Daniela Puiu, Alison G. Compton, Shannon Cowie, John Christodoulou, David R. Thorburn, Hongyu Zhao, Dan E. Arking, Shamil R. Sunyaev, Monkol Lek. Quantifying constraint in the human mitochondrial genome. bioRxiv 2023; doi: https://doi.org/10.1101/2022.12.16.520778