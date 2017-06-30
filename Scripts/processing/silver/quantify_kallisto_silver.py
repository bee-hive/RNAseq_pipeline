##################################################
#  quantify_kallisto_silver.py
#
#  $proj/Scripts/processing/silver/quantify_kallisto_silver.py
#
#  Create expression matrix for kallisto
#
#  Authors: Brian Jo
#
##################################################
import pandas as pd
from sys import argv

tables_dir = '/tigress/BEE/RNAseq/Data/Resources/gtex/tables/'
quant_dir = '/tigress/BEE/gtex/data/phenotype/expression/mapped_rna_seq_reads/silver/kallisto_hg38/'
out_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/silver/kallisto/'

sample_table = pd.read_csv(tables_dir + 'sample_table.txt', sep='\t')
subject_table = pd.read_csv(tables_dir + 'subject_table.txt', sep='\t')
tissue_table = pd.read_csv(tables_dir + 'tissue_table.txt', sep='\t')

# which subjects/samples have genotypes?
subjects_with_geno = list(subject_table[subject_table['genotype_avail'] == True]['submitted_subject_id_s'])
subjects_without_geno = list(subject_table[subject_table['genotype_avail'] == False]['submitted_subject_id_s'])

samples_with_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_with_geno)]
samples_without_geno = sample_table[sample_table['submitted_subject_id_s'].isin(subjects_without_geno)]
isoform_suffix = 'abundance.tsv'

# create table for each tissue
# for tissue in tissue_table['tissue_name']:
tissue = 'ovary'
print(tissue)
sample_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['Run_s'])
sample_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['Run_s'])
subject_list_with_geno = list(samples_with_geno[samples_with_geno['tissue_name'] == tissue]['submitted_subject_id_s'])
subject_list_without_geno = list(samples_without_geno[samples_without_geno['tissue_name'] == tissue]['submitted_subject_id_s'])

# gene_expr_mat_TPM = pd.DataFrame()
# gene_expr_mat_counts = pd.DataFrame()
isoform_expr_mat_TPM = pd.DataFrame()
isoform_expr_mat_counts = pd.DataFrame()
for i in range(len(sample_list_with_geno)):
    print(str(i))
    sample = sample_list_with_geno[i]
    subject = subject_list_with_geno[i]
    # gene_table = pd.read_csv(quant_dir + sample + '/' + sample + gene_suffix, sep='\t')
    isoform_table = pd.read_csv(quant_dir + sample + '/' + isoform_suffix, sep='\t')
    
    isoform_expr_TPM = isoform_table[['tpm']]
    isoform_expr_counts = isoform_table[['est_counts']]
    # name expression column by subject
    isoform_expr_TPM.columns = [subject]
    isoform_expr_counts.columns = [subject]
    # build entire expression matrix incrementally by column (subject)
    isoform_expr_mat_TPM = pd.concat([isoform_expr_mat_TPM, isoform_expr_TPM], axis=1)
    isoform_expr_mat_counts = pd.concat([isoform_expr_mat_counts, isoform_expr_counts], axis=1)

    # TODO: combine across genes

    # program may use a string like '.' to signify 0, in which case pandas will turn it into NA
    # expr_mat_TPM_uncertain = expr_mat_TPM_uncertain.dropna(axis=0)
    # expr_mat_TPM_certain = expr_mat_TPM_certain.dropna(axis=0)
    # expr_mat_counts_uncertain = expr_mat_counts_uncertain.dropna(axis=0)
    # expr_mat_counts_certain = expr_mat_counts_certain.dropna(axis=0)

# write files
isoform_expr_mat_TPM.index = isoform_table['target_id'].tolist()
isoform_expr_mat_counts.index = isoform_table['target_id'].tolist()

# Need to convert to gene table first and then Take out rows that don't have at least 10 nonzero values
n_zeros = (isoform_expr_mat_TPM == 0).sum(1)
rows_to_keep = n_zeros <= (isoform_expr_mat_TPM.shape[1] - 10)

isoform_expr_mat_TPM = isoform_expr_mat_TPM.loc[rows_to_keep]
isoform_expr_mat_counts = isoform_expr_mat_counts.loc[rows_to_keep]

# gene_expr_mat_TPM.to_csv(out_dir + 'genes/' + tissue + '_kallisto_TPM_with_genotype.txt', sep='\t')
# gene_expr_mat_counts.to_csv(out_dir + 'genes/' + tissue + '_kallisto_counts_with_genotype.txt', sep='\t')
isoform_expr_mat_TPM.to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_TPM_with_genotype.txt', sep='\t', float_format='%.4f')
isoform_expr_mat_counts.to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_counts_with_genotype.txt', sep='\t', float_format='%.4f')

    # If there are samples without genotype:
    # gene_expr_mat_TPM = pd.DataFrame()
    # gene_expr_mat_counts = pd.DataFrame()
    isoform_expr_mat_TPM = pd.DataFrame()
    isoform_expr_mat_counts = pd.DataFrame()
    for i in range(len(sample_list_without_geno)):
        print(str(i))
        sample = sample_list_without_geno[i]
        subject = subject_list_without_geno[i]
        # gene_table = pd.read_csv(quant_dir + sample + '/' + sample + gene_suffix, sep='\t')
        isoform_table = pd.read_csv(quant_dir + sample + '/' + isoform_suffix, sep='\t')
        
        isoform_expr_TPM = isoform_table[['tpm']]
        isoform_expr_counts = isoform_table[['est_counts']]
        # name expression column by subject
        isoform_expr_TPM.columns = [subject]
        isoform_expr_counts.columns = [subject]
        # build entire expression matrix incrementally by column (subject)
        isoform_expr_mat_TPM = pd.concat([isoform_expr_mat_TPM, isoform_expr_TPM], axis=1)
        isoform_expr_mat_counts = pd.concat([isoform_expr_mat_counts, isoform_expr_counts], axis=1)

        # TODO: combine across genes

        # program may use a string like '.' to signify 0, in which case pandas will turn it into NA
        # expr_mat_TPM_uncertain = expr_mat_TPM_uncertain.dropna(axis=0)
        # expr_mat_TPM_certain = expr_mat_TPM_certain.dropna(axis=0)
        # expr_mat_counts_uncertain = expr_mat_counts_uncertain.dropna(axis=0)
        # expr_mat_counts_certain = expr_mat_counts_certain.dropna(axis=0)

    # write files
    # gene_expr_mat_TPM.to_csv(out_dir + 'genes/' + tissue + '_kallisto_TPM_without_genotype.txt', sep='\t')
    # gene_expr_mat_counts.to_csv(out_dir + 'genes/' + tissue + '_kallisto_counts_without_genotype.txt', sep='\t')
    isoform_expr_mat_TPM.to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_TPM_without_genotype.txt', sep='\t')
    isoform_expr_mat_counts.to_csv(out_dir + 'isoforms/' + tissue + '_kallisto_counts_without_genotype.txt', sep='\t')

# Need to fix: 
# annotate the isoform/gene names on the first column
# remove rows that have no expression
# rounding error

