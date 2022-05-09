import argparse
import csv
import datetime
import numpy as np
import os
from typing import Dict, List, Tuple, Union


def initialize_sum_dict(identifier_list: List[str]):
	"""Generate a dictionary that will be used to sum the observed maximum heteroplasmy, mutation likelihood scores and
	count for a group of variants.

	:param identifier_list: list of the identifiers for each of the variant groups for which the values are being summed
	:return: a dictionary with a tuple key of the identifier, mutation group, value to sum, and region with variant,
	and value of 0 for each key
	"""
	# initialize dictionary so all values are 0
	dict = {}
	for identifier in identifier_list:
		for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to fit
			for value in ['obs_max_het', 'sum_LR', 'count']:  # the three values to sum
				for region in ['ref_exc_ori', 'ori']:  # the two mutational models that could be applied
					dict[identifier, mut_group, value, region] = 0
	return dict


def sum_obs_likelihood(
		mutation: str, identifier: str, region: str, observed: str, likelihood: str,
		dict: Dict[Tuple[str, str, str, str], Union[float, int]]):
	"""Sum the observed maximum heteroplasmy, mutation likelihood scores and count for a group of variants.

	:param mutation: the mutation type, in Ref>Alt format
	:param identifier: the identifier of the variant group for which the values are being summed
	:param region: should be either 'ref_exc_ori' or 'ori', to indicate which mutational model to apply
	:param observed: the value to use for the observed maximum heteroplasmy
	:param likelihood: the value to use for the mutation likelihood score
	:param dict: the name of the dictionary to append to
	:return: a dictionary with a tuple key of the identifier, mutation group, value to sum, and region with variant,
	and value corresponding to the value to sum in key
	"""
	# highly mutable G>A and T>C are fit separately
	mut_group = 'G>A_and_T>C' if (mutation == 'G>A' or mutation == 'T>C') else 'other'
	# note the dictionary should already be initialized (to start at 0 for all anticipated keys)
	dict[(identifier, mut_group, 'obs_max_het', region)] += float(observed)
	dict[(identifier, mut_group, 'sum_LR', region)] += float(likelihood)
	dict[(identifier, mut_group, 'count', region)] += 1
	return dict


def calibrate(input_file: str):
	"""Sum the observed maximum heteroplasmy in gnomAD, mutation likelihood scores and count for neutral variants,
	in the reference sequence excluding the OriB-OriH region.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	"""
	# exclude “artifact_prone_sites” in gnomAD positions - 301, 302, 310, 316, 3107, and 16182 (3107 already excluded)
	# variants at these sites were not called and therefore these positions removed from calculations
	artifact_prone_sites = [301, 302, 310, 316, 16182]
	
	# first, extract list of genes/loci and their length
	gene_length = {}
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		for gene in row["symbol"].split(','):  # split is to handle positions in two genes
			# mark control region separately to other non-coding
			gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
			gene = 'other_non-coding' if gene == '' else gene
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in artifact_prone_sites):
				if gene not in gene_length:
					gene_length[gene] = 1
				else:
					gene_length[gene] += 1
	
	# now initialize a dictionary that will be used to sum values for each gene/loci, and set all values to 0
	calibration = initialize_sum_dict(identifier_list=list(gene_length.keys()))
	
	# determine the phyloP threshold that will be used to identify neutral variants - this is the bottom decile
	phylop = []
	catch_list = []  # use to get to unique positions, since this file has three rows per positions (ie 3 possible SNVs)
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) not in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # second element [1] is bottom decile
	
	# now build dictionary
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for gene in row["symbol"].split(','):  # to handle variants within two genes
			if (int(row["POS"]) not in ori_region) and (int(row["POS"]) not in artifact_prone_sites):
				# criteria for neutral = haplogroup variant in phylotree or in bottom decile phyloP (threshold)
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# mark control region separately to other non-coding
					gene = 'control_region' if (int(row["POS"]) <= 576 or int(row["POS"]) >= 16024) else gene
					gene = 'other_non-coding' if gene == '' else gene
					# sum values for each gene/loci
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=gene, region='ref_exc_ori', dict=calibration,
						observed=row["gnomad_max_hl"], likelihood=row["Likelihood"])
	
	# write to file for plotting
	f = open('output_files/calibration/loci_obs_vs_scores.txt', "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to calibrate
		for gene in gene_length:
			f.write(
				mut_group + '\t' + gene + '\t' +
				str(calibration[(gene, mut_group, 'obs_max_het', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'sum_LR', 'ref_exc_ori')]) + '\t' +
				str(calibration[(gene, mut_group, 'count', 'ref_exc_ori')]) + '\t' +
				str(gene_length[gene] / 3) + '\n')  # need to divide by 3 as count is number positions x 3 (possible SNVs)


def calibrate_ori(input_file: str):
	"""Sum the observed maximum heteroplasmy in gnomAD, mutation likelihood scores and count for neutral variants,
	in the OriB-OriH region.

	:param input_file: annotated file with mutation likelihood scores and observed maximum heteroplasmy
	"""
	# first, segment the ori region into equally sized blocks of 70 bp, using approximate size of tRNA genes
	ori_blocks = {}
	block = 1
	for i in range(0, len(ori_region), 70):
		if len(ori_region[i:i + 70]) == 70:
			ori_blocks["block_" + str(block)] = ori_region[i:i + 70]
		else:  # to ensure all ori bases included, make the last block larger if needed
			ori_blocks["block_" + str(block - 1)] = ori_blocks["block_" + str(block - 1)] + \
													ori_region[i:i + ori_region[-1]]
		block += 1
	
	# now initialize a dictionary that will be used to sum values for each block, and set all values to 0
	calibration = initialize_sum_dict(identifier_list=list(ori_blocks.keys()))
	
	# determine the phyloP threshold that will be used to identify neutral variants - this is the bottom decile
	phylop = []
	catch_list = []  # use to get to unique positions, since this file has three rows per positions (ie 3 possible SNVs)
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		if int(row["POS"]) in ori_region:
			if row["POS"] not in catch_list:
				phylop.append(float(row["phyloP_score"]))
				catch_list.append(row["POS"])
	phylop_threshold = np.percentile(np.array(phylop), np.arange(0, 100, 10))[1]  # second element [1] is bottom decile
	
	# now build dictionary
	for row in csv.DictReader(open(input_file), delimiter='\t'):
		mutation = row["REF"] + '>' + row["ALT"]
		for n_block in ori_blocks:
			if int(row["POS"]) in ori_blocks[n_block]:  # if the variant is in the ori block in the loop
				# criteria for neutral = haplogroup variant in phylotree or in bottom decile phyloP (threshold)
				if int(row["in_phylotree"]) == 1 or (float(row["phyloP_score"]) < float(phylop_threshold)):
					# sum values for each block
					calibration = sum_obs_likelihood(
						mutation=mutation, identifier=n_block, region='ori', dict=calibration,
						observed=row["gnomad_max_hl"], likelihood=row["Likelihood"])
	
	# write to file for plotting
	f = open('output_files/calibration/loci_obs_vs_scores_ori.txt', "w")
	header = "mutation_group	symbol	obs_max_het	sum_likelihood	count	length"
	f.write(header + '\n')
	for mut_group in ['G>A_and_T>C', 'other']:  # the two mutation groups to calibrate
		for n_block in ori_blocks:
			f.write(
				mut_group + '\t' + n_block + '\t' +
				str(calibration[(n_block, mut_group, 'obs_max_het', 'ori')]) + '\t' +
				str(calibration[(n_block, mut_group, 'sum_LR', 'ori')]) + '\t' +
				str(calibration[(n_block, mut_group, 'count', 'ori')]) + '\t' +
				str(len(ori_blocks[n_block])) + '\n')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"-input", type=str, help="Annotated file with mutation likelihood scores and observed maximum heteroplasmy")
	args = parser.parse_args()
	
	# set default
	if args.input is None:
		args.input = 'output_files/mutation_likelihoods/mitochondrial_mutation_likelihoods_annotated.txt'
	
	for path in ['output_files/calibration']:
		if not os.path.exists(path):
			os.makedirs(path)
			print("Creating required directories")
	
	print(datetime.datetime.now(), "Preparing files for model calibration!" + '\n')
	
	# ori refers to OriB-OriH region with known difference in mutational signature, m.16197-191, across artificial break
	ori_region = list(range(16197, 16570)) + list(range(1, 191 + 1))
	
	calibrate(input_file=args.input)
	calibrate_ori(input_file=args.input)
	
	print(datetime.datetime.now(), "Script complete!" + '\n')
