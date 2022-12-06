## Overview:

`databases`: Files and their sources include:
- gnomAD dataset [from the gnomAD browser](https://gnomad.broadinstitute.org/downloads). 
- HelixMTdb dataset [from the Helix website](https://www.helix.com/pages/mitochondrial-variant-database).
- MITOMAP dataset [from MITOMAP](https://www.mitomap.org/MITOMAP/resources).
- Disease-associated variants [from MITOMAP](https://www.mitomap.org/MITOMAP/resources).
- Disease-associated variants [from ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/).
- Haplogroup-associated variants extracted [from PhyloTree](https://www.phylotree.org/).
- List of loci in the mitochondrial genome [from MITOMAP](https://www.mitomap.org/foswiki/bin/view/MITOMAP/GenomeLoci).

`input_denovo`: Mitochondrial de novo mutations ascertained from publication and in-house datasets. Details on their curation is provided in the Supplementary Information, see <placeholder>.

`insilicos`: 
- PhyloP conservation scores from alignment of 100 vertebrates, [from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.100way.phyloP100way/).
- APOGEE prediction scores [from MitImpact](https://mitimpact.css-mendel.it/).
- MitoTip prediction scores [from MITOMAP](https://www.mitomap.org/MITOMAP/MitoTipScores).
- HmtVar disease scores retrieved [from HmtVar](https://www.hmtvar.uniba.it/).

`other_annotations`:
- UniProt Knowledgebase protein annotations obtained [from the UniProt FTP site](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/). Note annotations for mtDNA-encoded proteins were extracted from all available bed files for Homo sapiens (available in UniProt directory *UP000005640_9606_beds*).
- Residues involved in complex I proton transfer curated from [PMID:32972993](https://pubmed.ncbi.nlm.nih.gov/32972993/).
- tRNA position numbers, source file obtained [from MitoTIP](https://github.com/sonneysa/MitoTIP/) and manually converted to a list of tRNA positions for each mtDNA coordinate.
- RNA modifications and tRNA domains, obtained as previously described [by Lake et al](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356).
- RNA base type extracted from manually curated data reported previously [by Lake et al](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356), using file `all_RNA_bases.tsv`.
- Bases involved in rRNA:rRNA bridges curated from [PMID:25838379](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4501431/)

`synthetic_vcf`: Files in this directory obtained as described [previously](https://github.com/broadinstitute/gnomad-mitochondria/tree/main/gnomad_mitochondria/manuscript_analyses).

`svg_input`: Files used to generate svg figures displaying the canonical secondary structure of a tRNA, or each rRNA, obtained as previously described [by Lake et al](https://academic.oup.com/bioinformatics/article/38/10/2967/6567356).
