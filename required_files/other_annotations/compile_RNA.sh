# generate a single file to be used to determine base type in RNA
head -1 /Users/nicole.lake/Dropbox/Yale/PROJECTS/MitoVisualize/2019_Summer_Student_Project/github_mito_SVGs/tsvs_final/MT-TA_with_pairs.tsv | tr -d '\r' | awk 'BEGIN{OFS="\t"} {print $0,"file"}' > required_files/other_annotations/all_RNA_bases.tsv

for file in /Users/nicole.lake/Dropbox/Yale/PROJECTS/MitoVisualize/2019_Summer_Student_Project/github_mito_SVGs/tsvs_final/*tsv;
do
	name=${file#*final/}
	echo $name
	tail -n +2 -q $file | tr -d '\r' | awk -v n=$name 'BEGIN{OFS="\t"} {print $0,n}' >> required_files/other_annotations/all_RNA_bases.tsv
done