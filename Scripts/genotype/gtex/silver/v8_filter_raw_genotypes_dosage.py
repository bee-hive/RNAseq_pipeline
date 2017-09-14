##################################################
#  v8_filter_raw_genotypes_dosage.py
#
#  $proj/Scripts/genotype/gtex/silver/v8_filter_raw_genotypes_dosage.py
# 
#  This script process the raw WGS vcf file and splits them into chromosomes, while filtering for MAF and indels
#
#  Author: Brian Jo
#
##################################################

# function for processing each line
def process_vcf_line(line):
	# We will apply the two following filters:
	# Filter based on MAF
	# Only take SNPs - filter out indels
	entry = line.strip().decode("utf-8").split('\t')
	processed_line = {}
	processed_line['chr'] = entry[0][3:]
	# don't process for MAF < 0.01
	processed_line['MAF'] = float(entry[7].split(';')[1][3:])
	if processed_line['MAF'] < 0.01:
		return processed_line
	processed_line['ID'] = entry[2]
	processed_line['SNP'] = True
	# don't process for indels
	if len(entry[3]) > 1 or len(entry[4]) > 1:
		processed_line['SNP'] = False
		return processed_line
	dosage_arr = [str(int(entry[x][0]) + int(entry[x][2])) if entry[x][0] != '.' else '-' for x in range(9, len(entry))]
	processed_line['dosage'] = dosage_arr
	return processed_line

import gzip
import os

vcf_raw = '/tigress/BEE/gtex/dbGaP-7716/57610/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz'
f = gzip.open(vcf_raw, 'rb')

while True:
	line = f.readline()
	# traverse until the suject IDs are found
	if line.strip().decode("utf-8").split('\t')[0] == '#CHROM':
		break

# This line contains the subject IDs
subject_IDs = line.strip().decode("utf-8").split('\t')[9:]
# len(subject_IDs) = 838
chrom_list = [str(x+1) for x in range(22)] + ['X', 'Y', 'M']

# We will make two files - for MAF 5% and 1%
out_dir = '/tigress/BEE/RNAseq/Data/Genotype/gtex/v8/allelic_dosage/'

current_chr = ''
# dummy files
out_f_01 = open(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_01.txt', 'w')
out_f_05 = open(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_05.txt', 'w')

for line in f:
	processed_line = process_vcf_line(line)
	# start new files:\
	if processed_line['chr'] not in chrom_list:
		continue
	if current_chr != processed_line['chr']:
		current_chr = processed_line['chr']
		print(current_chr)
		out_f_01.close()
		out_f_05.close()
		out_f_01 = open(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_01.txt', 'w')
		out_f_05 = open(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + current_chr + '_dosage_MAF_05.txt', 'w')
		# write header
		out_f_01.write('\t' + '\t'.join(subject_IDs) + '\n')
		out_f_05.write('\t' + '\t'.join(subject_IDs) + '\n')
	if processed_line['MAF'] >= 0.05 and processed_line['SNP']:
		out_f_05.write(processed_line['ID'] + '\t' + '\t'.join(processed_line['dosage']) + '\n')
	elif processed_line['MAF'] >= 0.01 and processed_line['SNP']:
		out_f_01.write(processed_line['ID'] + '\t' + '\t'.join(processed_line['dosage']) + '\n')

out_f_01.close()
out_f_05.close()

# remove dummy files
os.remove(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_dosage_MAF_01.txt')
os.remove(out_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr_dosage_MAF_01.txt')