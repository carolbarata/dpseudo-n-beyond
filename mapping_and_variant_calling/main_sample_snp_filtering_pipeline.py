##################################################
# Python3 script that filters SNPs from VCF
# extracted data from a single sample
##################################################

##### Argument parser #####
import argparse, os, re
parser = argparse.ArgumentParser(description="SNP filtering pipeline")

parser.add_argument("-b","--input_bwa_file",
required = True,
dest = "b",
type = str,
help = "Full path to input tab separated data file for bwa mem mapped reads")

parser.add_argument("-n","--input_novoalign_file",
required = True,
dest = "n",
type = str,
help = "Full path to input tab separated data file for novoalign mapped reads")

parser.add_argument("-t","--time_point",
required = True,
dest = "t",
type = int,
help = "Time point in your time series")

parser.add_argument("-s","--sample_size",
required = True,
dest = "s",
type = int,
help = "Pooled sequencing sample size")

parser.add_argument("-o","--output_file",
required = True,
dest = "o",
type = str,
help = "Full path to output file")

args = parser.parse_args()


#### Import required modules #####

import csv, math, os.path
import pandas as pd

def load_novoalign_file(novoalign_file):
	#### Read in tab sep data files ####
	#snp_file_bwa = open(bwa_file)
	#snp_data_bwa = csv.reader(snp_file_bwa, delimiter = '\t')
	
	snp_file_novoalign = open(novoalign_file)
	snp_data_novoalign = csv.reader(snp_file_novoalign, delimiter = '\t')
	
	##### Set up dictionary keys #####
	#headers_bwa = next(snp_data_bwa, None)
	#snp_table_bwa = {}
	#for variable in sorted(headers_bwa):
	#	snp_table_bwa[variable] = []
	
	headers_novoalign = next(snp_data_novoalign, None)
	snp_table_novoalign = {}
	for variable in sorted(headers_novoalign):
		snp_table_novoalign[variable] = []
	
	##### Fill in dictionary values #####
	#for snp in snp_data_bwa:
	#	for variable, entry in zip(headers_bwa, snp):
	#		snp_table_bwa[variable].append(entry)
	
	for snp in snp_data_novoalign:
		for variable, entry in zip(headers_novoalign, snp):
			snp_table_novoalign[variable].append(entry)
	#print('Loading done') if len(snp_table_bwa['POS']) > 1 and len(snp_table_novoalign['POS']) > 1 else print('Loading failed')
	return snp_table_novoalign


def filter_indexes(snp_entry_bwa, snp_table_novoalign, time_point, sample_size):
	##### Compute strand bias #####
	#snp_table_bwa['ADF'] = [str(x) for x in snp_table_bwa['ADF']]
	#snp_table_novoalign['ADF'] = [str(x) for x in snp_table_novoalign['ADF']]
	#tot_adf_bwa = [x.split(",") for x in snp_table_bwa['ADF']]
	#tot_adf_novoalign = [x.split(",") for x in snp_table_novoalign['ADF']]
	tot_adf_bwa = snp_entry_bwa['ADF'][0].split(",") if ',' in snp_entry_bwa['ADF'][0] else snp_entry_bwa['ADF'][0]
	#tot_adf_novoalign = snp_table_novoalign['ADF'].split(",") if ',' in snp_table_novoalign['ADF'] else snp_table_novoalign['ADF']

	#### Filter based on multi-/biallelic status ####
	#multiallele_index_bwa = [index for index, snp in enumerate(tot_adf_bwa) if len(snp) == 2]
	multiallele_index_bwa = True if len(tot_adf_bwa) == 2 else False
	
	#tot_adf_bwa = [int(float(x[0])) + int(float(x[1])) if len(x) >= 2 else int(float(x)) for x in tot_adf_bwa]
	tot_adf_bwa = int(float(tot_adf_bwa[0])) + int(float(tot_adf_bwa[1])) if len(tot_adf_bwa) >= 2 else int(float(tot_adf_bwa[0]))
	#tot_adf_bwa = [x / int(float(y)) for x, y in zip(tot_adf_bwa, snp_table_bwa['DP']) if int(float(y)) > 0]
	tot_adf_bwa = tot_adf_bwa / int(float(snp_entry_bwa['DP'][0])) if int(float(snp_entry_bwa['DP'][0])) > 0 else 0
	#stn_bias_bwa = [abs(x - 0.5) for x in tot_adf_bwa]
	stn_bias_bwa = abs(tot_adf_bwa - 0.5)
	snp_entry_bwa['BIAS'] = stn_bias_bwa
	
	#tot_adf_novoalign = [int(float(x[0])) + int(float(x[1])) if len(x) >= 2 else int(float(x)) for x in tot_adf_novoalign]
	#tot_adf_novoalign = int(float(tot_adf_novoalign[0])) + int(float(tot_adf_novoalign[1])) if len(tot_adf_novoalign) >= 2 else int(float(tot_adf_novoalign))
	#tot_adf_novoalign = [x / int(float(y)) for x, y in zip(tot_adf_novoalign, snp_table_novoalign['DP']) if int(float(y)) != 0]
	#tot_adf_novoalign = tot_adf_novoalign / int(float(snp_table_novoalign['DP'])) if int(float(snp_table_novoalign['DP'])) > 0 else 0
	#stn_bias_novoalign = [abs(x - 0.5) for x in tot_adf_novoalign]
	#stn_bias_novoalign = abs(tot_adf_novoalign - 0.5)
	#snp_table_novoalign['BIAS'] = stn_bias_novoalign
	##### Filter to be applied on first time point samples #####
	if time_point == 1:
		#### Filter based on genotype ####
		#poly_gt_index_bwa = [index for index, genotype in enumerate(snp_table_bwa['GT'])  if genotype != '1/1']
		poly_gt_index_bwa = True if snp_entry_bwa['GT'][0] != '1/1' else False
	
	##### Compute ref/alt allele frequencies #####
	#snp_table_bwa['AD'] = [str(x) for x in snp_table_bwa['AD']]
	#ref_freq_bwa = [int(float(x.split(",")[0])) / int(float(y)) if ',' in x and int(float(y)) != 0 else int(float(x)) / int(float(y)) for x, y in zip(snp_table_bwa['AD'], snp_table_bwa['DP'])]
	ref_freq_bwa = int(float(snp_entry_bwa['AD'][0].split(",")[0])) / int(float(snp_entry_bwa['DP'][0])) if ',' in snp_entry_bwa['AD'][0] and int(float(snp_entry_bwa['DP'][0])) != 0 else int(float(snp_entry_bwa['AD'][0])) / int(float(snp_entry_bwa['DP'][0]))
	snp_entry_bwa['REF_FQ'] = ref_freq_bwa
	#alt_freq_bwa = [int(float(x.split(",")[1])) / int(float(y)) if ',' in x and int(float(y)) != 0 else 0 for x, y in zip(snp_table_bwa['AD'], snp_table_bwa['DP'])]
	alt_freq_bwa = int(float(snp_entry_bwa['AD'][0].split(",")[1])) / int(float(snp_entry_bwa['DP'][0])) if ',' in snp_entry_bwa['AD'][0] and int(float(snp_entry_bwa['DP'][0])) != 0 else int(float(snp_entry_bwa['AD'][0])) / int(float(snp_entry_bwa['DP'][0]))
	#snp_table_novoalign['AD'] = [str(x) for x in snp_table_novoalign['AD']]
	#ref_freq_novoalign = [int(float(x.split(",")[0])) / int(float(y)) if ',' in x and int(float(y)) != 0 else 0 for x, y in zip(snp_table_novoalign['AD'], snp_table_novoalign['DP'])]
	#ref_freq_novoalign = int(float(snp_table_novoalign['AD'].split(",")[0])) / int(float(snp_table_novoalign['DP'])) if ',' in snp_table_novoalign['AD'] and int(float(snp_table_novoalign['DP'])) != 0 else int(float(snp_table_novoalign['AD'])) / int(float(snp_table_novoalign['DP']))
	#snp_table_novoalign['REF_FQ'] = ref_freq_novoalign
	
	##### Filter based on min and max allele frequency #####
	min_af = 1 / sample_size # 1 / n where n is the sample size
	max_af = 1 - min_af
	#ref_af_filter_index_bwa = [index for index, allele_freq in enumerate(ref_freq_bwa) if allele_freq > min_af and allele_freq < max_af]
	ref_af_filter_index_bwa = True if ref_freq_bwa > min_af and ref_freq_bwa < max_af else False
	
	##### Filter based on multi-/biallelic status #####
	#sum_freq_bwa = [x + y for x, y in zip(ref_freq_bwa, alt_freq_bwa)]
	sum_freq_bwa = ref_freq_bwa + alt_freq_bwa
	#sum_freq_index_bwa = [index for index, total in enumerate(sum_freq_bwa) if total == 1]
	sum_freq_index_bwa = True if sum_freq_bwa == 1 else False
	
	##### Filter variants based on quality criteria #####
	if time_point == 1:
		#all_filter_index_bwa = list(set(multiallele_index_bwa).intersection(set(poly_gt_index_bwa), set(ref_af_filter_index_bwa), set(sum_freq_index_bwa)))
		all_filter_bwa = multiallele_index_bwa and poly_gt_index_bwa and ref_af_filter_index_bwa and sum_freq_index_bwa
	else:
		#all_filter_index_bwa = list(set(multiallele_index_bwa).intersection(set(ref_af_filter_index_bwa), set(sum_freq_index_bwa)))
		all_filter_bwa = multiallele_index_bwa and ref_af_filter_index_bwa and sum_freq_index_bwa
	
	return all_filter_bwa


def two_mapper_filter(snp_entry_bwa, snp_table_novoalign, filtered_index):
		
	##### Intersect bwa mem and novoalign mapped read files #####
	#common_snps = list(set(snp_data_bwa['POS']).intersection(set(snp_data_novoalign['POS'])))
	common_snp = snp_entry_bwa['POS'] in snp_table_novoalign['POS']
	#common_snp_index = [snp_data_bwa['POS'].index(snp) for snp in common_snps]
	#keeper_index = list(set(all_filter_index_bwa).intersection(set(common_snp_index)))
	
	#for key, values in snp_data_bwa.items():
	#	snp_data_bwa[key] = [i for j, i in enumerate(values) if j in keeper_index]
	
	#for index, snp in enumerate(snp_table_bwa['POS']):
	#	novoalign_snp_index = snp_table_novoalign['POS'].index(snp)
	#	if int(float(snp_data_bwa['DP'][index])) < int(float(snp_table_novoalign['DP'][novoalign_snp_index])):
	#		for key in snp_data_bwa.keys():
	#			snp_data_bwa[key][index] = snp_table_novoalign[key][novoalign_snp_index]
	if common_snp:
		common_snp_index = snp_table_novoalign['POS'].index(snp_entry_bwa['POS'])
		if int(float(snp_entry_bwa['DP'])) < int(float(snp_table_novoalign['DP'][common_snp_index])) and snp_entry_bwa['REF'] == snp_table_novoalign['REF'][common_snp_index] and snp_entry_bwa['ALT'] == snp_table_novoalign['ALT'][common_snp_index]:
			for key in snp_entry_bwa.keys():
				snp_entry_bwa[key] = snp_table_novoalign[key][common_snp_index]
	
	return snp_entry_bwa


def filter_snp_by_snp(snp_entry_bwa, snp_table_novoalign, time_point, sample_size, output_file):
	
	##### First round of filtering #####
	filtered_index = filter_indexes(snp_entry_bwa, snp_table_novoalign, time_point, sample_size)
	
	if filtered_index:
		##### Second round of filtering #####
		filtered_snp_entry = two_mapper_filter(snp_entry_bwa, snp_table_novoalign, filtered_index)
		pd_snp_entry = pd.DataFrame(filtered_snp_entry, index = [0]) 
		##### Save SNP data to existing csv #####
		if os.path.exists(output_file):
			pd_snp_entry.to_csv(output_file, sep = '\t', mode = 'a', header = False, index = False)
		##### Save SNP data to new csv #####
		else:
			pd_snp_entry.to_csv(output_file, sep = '\t', header = True, index = False)


def parse_novoalign_table(novoalign_file):
	#### Read in tab sep data files ####
	#snp_file_bwa = open(bwa_file)
	#snp_data_bwa = csv.reader(snp_file_bwa, delimiter = '\t')
	
	snp_file_novoalign = open(novoalign_file)
	snp_data_novoalign = csv.reader(snp_file_novoalign, delimiter = '\t')
	
	##### Set up dictionary keys #####
	#headers_bwa = next(snp_data_bwa, None)
	#snp_table_bwa = {}
	#for variable in sorted(headers_bwa):
	#       snp_table_bwa[variable] = []
	
	headers_novoalign = next(snp_data_novoalign, None)
	snp_table_novoalign = {}
	for variable in sorted(headers_novoalign):
		snp_table_novoalign[variable] = []
	
	##### Fill in dictionary values #####
	#for snp in snp_data_bwa:
	#       for variable, entry in zip(headers_bwa, snp):
	#               snp_table_bwa[variable].append(entry)
	
	for snp in snp_data_novoalign:
		for variable, entry in zip(headers_novoalign, snp):
			snp_table_novoalign[variable].append(entry)
	#print('Loading done') if len(snp_table_bwa['POS']) > 1 and len(snp_table_novoalign['POS']) > 1 else print('Loading failed')
	##### Compute strand bias #####
	snp_table_novoalign['ADF'] = [str(x) for x in snp_table_novoalign['ADF']]
	tot_adf_novoalign = [x.split(",") for x in snp_table_novoalign['ADF']]
	tot_adf_novoalign = [int(float(x[0])) + int(float(x[1])) if len(x) >= 2 else int(float(x[0])) for x in tot_adf_novoalign]
	tot_adf_novoalign = [x / int(float(y)) for x, y in zip(tot_adf_novoalign, snp_table_novoalign['DP']) if int(float(y)) != 0]
	stn_bias_novoalign = [abs(x - 0.5) for x in tot_adf_novoalign]
	snp_table_novoalign['BIAS'] = stn_bias_novoalign	
	##### Compute ref/alt allele frequencies #####
	snp_table_novoalign['AD'] = [str(x) for x in snp_table_novoalign['AD']]
	ref_freq_novoalign = [int(float(x.split(",")[0])) / int(float(y)) if ',' in x and int(float(y)) != 0 else 0 for x, y in zip(snp_table_novoalign['AD'], snp_table_novoalign['DP'])]
	snp_table_novoalign['REF_FQ'] = ref_freq_novoalign
		
	return snp_table_novoalign


def read_n_parse(bwa_file, snp_table_novoalign, time_point, sample_size, output_file):
	##### Open bwa data file #####
	with open(bwa_file) as snp_bwa_file:
		snp_count = 0
		##### Read one snp data line at a time #####
		for line in snp_bwa_file:
			snp_line_bwa = line.rstrip('\n').split('\t')
			if snp_count == 0:
				headers_bwa = snp_line_bwa
			else:
				print('Processing SNP', snp_count)
				snp_entry_bwa = {}
				for variable, entry in zip(headers_bwa, snp_line_bwa):
					snp_entry_bwa[variable] = []
					snp_entry_bwa[variable].append(entry)
				##### Parse snp data #####
				filter_snp_by_snp(snp_entry_bwa = snp_entry_bwa, snp_table_novoalign = snp_table_novoalign, time_point = time_point, sample_size = sample_size, output_file = output_file)
			snp_count += 1	


##### Parse novoalign data #####

snp_table_novoalign_psd = parse_novoalign_table(novoalign_file = args.n)

##### Read in and parse the data #####

read_n_parse(bwa_file = args.b, snp_table_novoalign = snp_table_novoalign_psd, time_point = args.t, sample_size = args.s, output_file = args.o)
	

