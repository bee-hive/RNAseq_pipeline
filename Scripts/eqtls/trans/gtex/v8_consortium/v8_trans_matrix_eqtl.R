##################################################
#  v8_trans_matrix_eqtl.R
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v8_trans_matrix_eqtl.R
# 
#  This version performs the all-by-all trans-mapping analysis for v8 consortium.
#
#  Author: Brian Jo
#
##################################################

#!/usr/bin/Rscript

# Manually import PATH
# .libPaths("/home/bj5/R/x86_64-redhat-linux-gnu-library/3.1")

args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# Example
args = c(1:7)
args[1] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/Whole_Blood.v8.normalized_expression.bed.gz'
args[2] = '1'
args[3] = '2'
args[4] = '100'
args[5] = '2.5e8'
args[6] = 'Whole_Blood'
args[7] = '/tigress/BEE/RNAseq/Data/Expression/gtex/hg38/GTEx_Analysis_v8_eQTL_expression_matrices/GTEx_Analysis_v8_eQTL_covariates/'
args[8] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/all-by-all-PEER-increments/skinnotsunexposedsuprapubic/skinnotsunexposedsuprapubic_nonverlapping_certain_autosomes_normalized_MatrixEQTL'

expression_file_location = args[1]
# changed feature - geno_option is always continuous by default
part_number = as.numeric(args[2])
chr_number = args[3]
partition_size = as.numeric(args[4])
cis_dist = as.numeric(args[5])
tissue_name = args[6]
cov_dir = args[7]
out_file = args[8]

# Read in the expression files and the gene positions
header = readLines(gzfile(expression_file_location), n = 1)
header = strsplit(header, '\t')[[1]]
expression_matrix = read.csv(gzfile(expression_file_location, 'r'), sep = '\t', stringsAsFactors = FALSE)

colnames(expression_matrix) = header
rownames(expression_matrix) = expression_matrix$gene_id

gene_positions = expression_matrix[,c(4,1,2,3)]
colnames(gene_positions) = c('gene_id', 'chr', 'start', 'end')
expression_matrix = expression_matrix[,c(5:ncol(expression_matrix))]

# Read in the genotype positions
genotype_file_name = paste0(proj_dir, '/Data/Genotype/gtex/v8/allelic_dosage/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr', chr_number, '_dosage_MAF_05.txt')
genotype_matrix_master = read.table(genotype_file_name, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
# Get the appropriate partition
num_parts = ceiling(nrow(genotype_matrix_master) / partition_size)
num_inds = ceiling(nrow(genotype_matrix_master) / num_parts)
if (part_number == num_parts) {
  genotype_matrix = genotype_matrix_master[c((((part_number-1)*partition_size)+1):nrow(genotype_matrix_master)),]
} else {
  genotype_matrix = genotype_matrix_master[c((((part_number-1)*partition_size)+1):(part_number*partition_size)),]
}

# There should be 838 indivs for genotype
rownames(genotype_matrix) = genotype_matrix$X
genotype_matrix = genotype_matrix[,c(2:ncol(genotype_matrix))]
colnames(genotype_matrix) = as.character(sapply(colnames(genotype_matrix), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))

# Make sure the columns are the same
expression_matrix = expression_matrix[,(colnames(expression_matrix) %in% colnames(genotype_matrix))]
genotype_matrix = genotype_matrix[,sapply(colnames(expression_matrix), function(x) {match(x, colnames(genotype_matrix))})]

# Fix the data type
genotype_matrix_temp = data.frame(lapply(genotype_matrix,as.numeric))
rownames(genotype_matrix_temp) = rownames(genotype_matrix)
genotype_matrix = genotype_matrix_temp

# Get the SNP positions
snp_positions = data.frame(ID = rownames(genotype_matrix), chr = sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][1]}), pos = as.numeric(sapply(rownames(genotype_matrix), function(x) {strsplit(x, '_')[[1]][2]})))

# Load in the covariates
suffix = '.v8.covariates.txt'
covars = read.csv(paste0(cov_dir, tissue_name, suffix), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
colnames(covars) = as.character(sapply(colnames(covars), function(x) {paste(strsplit(x, '\\.')[[1]][1], strsplit(x, '\\.')[[1]][2], sep = '-')}))
rownames(covars) = covars$ID
covars = covars[,colnames(expression_matrix)]

# Import MatrixEQTL
library(MatrixEQTL)
source(paste0(proj_dir, '/Scripts/eqtls/MatrixEQTL_wrapper.R'))

# Only for all-by-all runs: save the p-values that are over 1e-5, for which there are typically about 1e6 to 1e7 values
pvOutputThreshold = 1e-5
useModel = modelLINEAR;

me_trans_list = list()
MAF_list = list()
total_ntests = 0
hist_counts = rep(0,100)

for (i in c(1:nrow(genotype_matrix))) {
  print(i)
  inds = !is.na(genotype_matrix[i,])
  MAF_list[[i]] = sum(as.numeric(genotype_matrix[i,inds])) / (length(inds)*2)
  MAF_list[[i]] = min(MAF_list[[i]], 1 - MAF_list[[i]])

  snps = SlicedData$new()
  snps$CreateFromMatrix(as.matrix(genotype_matrix[i,inds]));
  cvrt = SlicedData$new()
  cvrt$CreateFromMatrix(as.matrix(covars[,inds]));
  gene = SlicedData$new()
  gene$CreateFromMatrix(as.matrix(expression_matrix[,inds]));

  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = NULL,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = pvOutputThreshold,
    snpspos = snp_positions[i,],
    genepos = gene_positions,
    cisDist = cis_dist,
    pvalue.hist = 100,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

  total_ntests = total_ntests + me$trans$ntests
  hist_counts = hist_counts + me$trans$hist.counts
  if (me$trans$neqtls > 0) {
    me_trans_list[[i]] = me$trans$eqtls
  }

}

# Check for existence of summary file to see if job completed
write.table(summary, file=paste0(out_file, '_part', sprintf("%03d", part_number), '_summary.txt'), row.names = FALSE, quote = FALSE, sep = '\t')
