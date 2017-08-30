##################################################
#  v8_normalization_methods_wrapper.py
#
#  $proj/Scripts/processing/covariates_exploratory/v8_normalization_methods_wrapper.py
#
#  Trying out various normalization methods/factor analysis in v8 expression data, without the known covariates
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd
import subprocess as sp
from sys import argv
import os

proj_dir = os.environ['proj']

tables_dir = proj_dir + '/Data/Resources/gtex/tables/v8/'
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

master_script = proj_dir + "/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods.sh"
master_handle = open(master_script, 'w')
master_handle.write("#!/bin/bash\n\n")

tissue_table.index = tissue_table['tissue_name']

n_factors_array = [10, 20, 30, 40, 50]

methods = argv[1:]

# Script for SFAmix
# We will try out gene qn and count matrices
if 'sfamix' in methods:
	out_dir = proj_dir + '/Output/processing/exploratory/v8/sfamix/no_cov/'
	exp_dir = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/'
	for tissue in tissue_table.index:
		if tissue_table.loc[tissue]['num_samples'] > 10:
			sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_sfamix_' + tissue + '.sh'
			sbatchhandle=open(sbatchfile, 'w')
			header=r"""#!/bin/bash
#SBATCH -J %s_sfamix      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/%s_sfamix
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue)
			sbatchhandle.write(header)
			# Is the inverted matrix in scratch?
			exp_mat = 'v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
			if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
				# Script for inverting the matrix
				sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/quantile_norm/ ' + exp_mat + '\n')
			# Is the inverted matrix in scratch?
			exp_mat = 'v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
			if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
				# Script for inverting the matrix
				sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/log_transform/ ' + exp_mat + '\n')
			for n_factors in n_factors_array:
				# Create directory
				directory = out_dir + 'quantile_norm/sample_norm/' + str(n_factors)
				if not os.path.exists(directory):
					os.makedirs(directory)
				directory = directory + '/' + tissue
				if not os.path.exists(directory):
					os.makedirs(directory)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				# Create directory
				directory = out_dir + 'log_transform/' + str(n_factors)
				if not os.path.exists(directory):
					os.makedirs(directory)
				directory = directory + '/' + tissue
				if not os.path.exists(directory):
					os.makedirs(directory)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
			sbatchhandle.close()
			master_handle.write("sbatch " + sbatchfile + '\n')

# Script for PEER factors
if 'peer' in methods:
	out_dir = proj_dir + '/Output/processing/exploratory/v8/sfamix/no_cov/'
	exp_dir = '/tigress/BEE/gtex/data/phenotype/expression/expression_matrices/v8/'
	for tissue in tissue_table.index:
		if tissue_table.loc[tissue]['num_samples'] > 10:
			sbatchfile = proj_dir + '/Scripts/processing/covariates_exploratory/batch/v8_normalization_methods_sfamix_' + tissue + '.sh'
			sbatchhandle=open(sbatchfile, 'w')
			header=r"""#!/bin/bash
#SBATCH -J %s_sfamix      # job name
#SBATCH -o %s/Scripts/processing/covariates_exploratory/joblogs/%s_sfamix
#SBATCH --get-user-env
#SBATCH --mem=30000           # 30 GB requested
#SBATCH -N 1                  # default
#SBATCH --ntasks-per-node=1   # 1 is default
#SBATCH -t 24:00:00           # walltime is an argument

# set default permissions on files created to be group rwx
umask 002

"""%(tissue, proj_dir, tissue)
			sbatchhandle.write(header)
			# Is the inverted matrix in scratch?
			exp_mat = 'v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
			if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
				# Script for inverting the matrix
				sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/quantile_norm/ ' + exp_mat + '\n')
			# Is the inverted matrix in scratch?
			exp_mat = 'v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
			if not os.path.exists('/scratch/gpfs/bj5/' + exp_mat):
				# Script for inverting the matrix
				sbatchhandle.write('Rscript ' + proj_dir + '/Scripts/processing/covariates_exploratory/v8_sfamix_invert_matrix.R ' + exp_dir + '/log_transform/ ' + exp_mat + '\n')
			for n_factors in n_factors_array:
				# Create directory
				directory = out_dir + 'quantile_norm/sample_norm/' + str(n_factors)
				if not os.path.exists(directory):
					os.makedirs(directory)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_tpm_' + tissue + '_sample_norm.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
				# Create directory
				directory = out_dir + 'log_transform/' + str(n_factors)
				if not os.path.exists(directory):
					os.makedirs(directory)
				exp_mat = '/scratch/gpfs/bj5/v8_RSEMv1.3.0_gene_count_' + tissue + '_log_transform.txt'
				param_array = ['/tigress/BEE/bin/SFAmix_code_documentation/SFAmix', '--nf', str(n_factors), '--y', exp_mat, '--out', directory, '--sep', 'tab']
				sbatchhandle.write(' '.join(param_array) + '\n')
			sbatchhandle.close()
			master_handle.write("sbatch " + sbatchfile + '\n')


print("sh " + master_script)
master_handle.close()
