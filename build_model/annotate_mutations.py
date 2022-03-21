import argparse
from compile_denovo import rcrs_pos_to_ref
import csv
import datetime


# note - add each annotation as needed downstream, so more to be added


def rcrs_pos_to_trinucleotide():
	"""Generate dictionary linking each position to its reference trinucleotide in the rCRS.

	:return: dictionary where the key is the position in rCRS, and the value is its reference trinucleotide
	"""
	# first, generate dictionary to convert coordinate to reference nucleotide
	rcrs_pos2ref = rcrs_pos_to_ref()
	# now, generate dictionary of coordinate to trinucleotide
	dict = {}
	for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'), delimiter='\t'):
		pos = int(row["POS"])
		ref = row["REF"]

		if pos == 16569:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(1)]  # dealing with circular genome
		elif pos == 1:
			trinucleotide = rcrs_pos2ref[str(16569)] + ref + rcrs_pos2ref[str(pos + 1)]  # dealing with circular genome
		elif pos == 3107:  # to handle the 'N'
			continue  # ie skip
		elif pos == 3106:  # to handle the 'N'
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(3108)]
		elif pos == 3108:  # to handle the 'N'
			trinucleotide = rcrs_pos2ref[str(3106)] + ref + rcrs_pos2ref[str(pos + 1)]
		else:
			trinucleotide = rcrs_pos2ref[str(pos - 1)] + ref + rcrs_pos2ref[str(pos + 1)]
		dict[str(pos)] = trinucleotide
	return dict


def gnomad_annotate():
	"""Generate dictionary with the maximum observed heteroplasmy (max_hl) for each variant in gnomAD.

	:return: a dictionary, where the key is a tuple of the position, ref and alt, and the value is the max_hl
	"""
	dict = {}
	for row in csv.DictReader(
			open('required_files/databases/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv'), delimiter='\t'):
		if row["filters"] == "PASS":
			dict[(row["position"], row["ref"], row["alt"])] = row["max_observed_heteroplasmy"]
	return dict


def phylop_annotate():
	"""Generate dictionary with the PhyloP scores for conservation, from 100 vertebrates.

	:return: a dictionary, where the key is the position, and the value is the PhyloP conservation score
	"""
	dict = {}
	pos = 0  # to handle a header row
	for row in open('required_files/insilicos/chrM.phyloP100way.wigFix'):
		dict[pos] = row.replace('\n', '')
		pos += 1
	return dict


def vep_annotate():
	"""Create a dictionary of the VEP annotations for every possible single nucleotide variant in the mtDNA.

	:return: dictionary where tuple of the variant and value identifier is key, and value is list of annotations
	"""
	vep = {}
	# use vcf where variants in two genes are split across two rows, for easy parsing
	for row in csv.DictReader(
			open("required_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf"), delimiter="\t"):
		# gene or locus
		if (row["REF"], row["POS"], row["ALT"], "symbol") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")].append(row["SYMBOL"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "symbol")] = [row["SYMBOL"]]
		# consequences
		if (row["REF"], row["POS"], row["ALT"], "consequence") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")].append(row["Consequence"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "consequence")] = [row["Consequence"]]
		# amino acids
		if (row["REF"], row["POS"], row["ALT"], "aa") in vep:  # i.e. variant is in two genes
			vep[(row["REF"], row["POS"], row["ALT"], "aa")].append(row["Amino_acids"])
		else:
			vep[(row["REF"], row["POS"], row["ALT"], "aa")] = [row["Amino_acids"]]
	return vep


def tRNA_positions():
	"""Create a dictionary of the tRNA positions (ranging from 1-73) encoded by mtDNA.

	:return: dictionary where position in the mtDNA is key and the value is the tRNA position
	"""
	dict = {}
	for row in csv.DictReader(open('required_files/other_annotations/tRNA_positions.txt'), delimiter='\t'):
		if row["m.pos"] in dict:  # if the mtDNA position is in two tRNA genes
			dict[row["m.pos"]].append(row["tRNA_position"])
		else:
			dict[row["m.pos"]] = [row["tRNA_position"]]
	return dict


def annotate(input_file: str):
	"""Annotate the file with all possible mitochondrial mutations and their likelihood scores.

	:param input_file: the file with mutation likelihood scores, output of composite_likelihood_mito.py
	"""
	f = open('%s_annotated.txt' % input_file.split('.txt')[0], "w")
	header = "POS	REF	ALT	Likelihood	trinucleotide	symbol	consequence	amino_acids	gnomad_max_hl	in_phylotree	phyloP_score	tRNA_position"
	f.write(header + '\n')

	# generate required dictionaries
	rcrs_pos2trinuc = rcrs_pos_to_trinucleotide()
	gnomad = gnomad_annotate()
	phylop = phylop_annotate()
	vep = vep_annotate()
	tRNA_position = tRNA_positions()
	
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		variant = row["REF"] + row["POS"] + row["ALT"]
		# annotate haplogroup variants, use '\n' in search to ensure unique matches
		in_phylo = 1 if ("\n" + variant + "\n") in open('required_files/databases/phylotree_variants.txt').read() else 0
		max_hl = gnomad[(row["POS"], row["REF"], row["ALT"])] if (row["POS"], row["REF"], row["ALT"]) in gnomad else 0
		tRNA_pos = tRNA_position[row["POS"]] if row["POS"] in tRNA_position else ''

		f.write(
			row["POS"] + '\t' + row["REF"] + '\t' + row["ALT"] + '\t' +
			row["Likelihood"] + '\t' +
			rcrs_pos2trinuc[row["POS"]] + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "symbol")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "consequence")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(vep[(row["REF"], row["POS"], row["ALT"], "aa")]).strip('[]').replace("'", "").replace(" ", "") + '\t' +
			str(max_hl) + '\t' +
			str(in_phylo) + '\t' +
			phylop[int(row["POS"])] + '\t' +
			str(tRNA_pos).strip('[]').replace("'", "").replace(" ", "") + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-input", type=str, help="File with mutation likelihood scores")
	args = parser.parse_args()

	# set default
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mitochondrial_mutation_likelihoods.txt'

	print(datetime.datetime.now(), "Annotating mitochondrial mutation likelihoods!" + '\n')

	annotate(input_file=args.input)

	print(datetime.datetime.now(), "Script complete!" + '\n')
