##################################################
#  v8_sfamix_incremental.R
#
#  $proj/Scripts/processing/covariates_exploratory/v8_sfamix_incremental.R
#
#  SFAmix takes in the inverted matrix - this script inverts the input matrix
#
#  Authors: Brian Jo
#
##################################################

library(SFAmix)
args <-commandArgs(TRUE)

# exp_dir = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/'
# output_dir = '/tigress/BEE/RNAseq/Output/processing/exploratory/v8/sfamix/GTEx_Analysis_v8/Adipose_Subcutaneous/'
# tissue = 'Adipose_Subcutaneous'
# n_factors = 500
# n_itr = 100

exp_dir = args[1]
output_dir = args[2]
tissue = args[3]
n_factors = args[4]
n_itr = args[5]

suffix = '.v8.normalized_expression.bed.gz'
output_file_suffix = paste0('_v8_', n_itr, 'sfamix_output.RData')
expression_file_location = paste0(exp_dir, tissue, suffix)

header = readLines(gzfile(expression_file_location), n = 1)
header = strsplit(header, '\t')[[1]]
expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = FALSE)

colnames(expression_matrix) = header
rownames(expression_matrix) = expression_matrix$gene_id
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

sfamix_result = SFAmixR(y = t(expression_matrix), nf = as.numeric(n_factors), itr = as.numeric(n_itr))

save(sfamix_result, file = paste0(output_dir, tissue, output_file_suffix))
