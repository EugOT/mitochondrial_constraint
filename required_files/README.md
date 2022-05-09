## Overview:

`databases`: Files include:
- gnomAD data in tsv format, obtained from the [gnomAD browser](https://gnomad.broadinstitute.org/downloads). 
- List of loci in mitochondrial genome [from MITOMAP](https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci).
- List of haplogroup variants from [PhyloTree](https://www.phylotree.org/), extracted from the htm file using in-house scripts.

`input_denovo`: Mitochondrial de novo mutations ascertained from publication and in-house datasets. Details on their curation is provided in the Supplementary Information, see <placeholder>.

`insilicos`: 
- PhyloP conservation scores from alignment of 100 vertebrates, obtained [from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way/).

`other_annotations`:
- tRNA positions, source file obtained [from MitoTIP](https://github.com/sonneysa/MitoTIP/) and manually converted to list of tRNA positions for each mtDNA coordinate.

`synthetic_vcf`: Files in this directory obtained as described [previously](https://github.com/broadinstitute/gnomad-mitochondria/tree/main/gnomad_mitochondria/manuscript_analyses).
