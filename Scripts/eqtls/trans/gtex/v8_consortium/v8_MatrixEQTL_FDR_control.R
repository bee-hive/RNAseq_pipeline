##################################################
#  v8_MatrixEQTL_FDR_control.R
#
#  $proj/Scripts/eqtls/trans/gtex/v8_consortium/v8_MatrixEQTL_FDR_control.R
# 
#  This version is the most up-to-date version for trans- pipeline.
#
#  Author: Brian Jo
#
##################################################

# Collect arguments
args <-commandArgs(TRUE)
proj_dir = Sys.getenv('proj')

# FDR Control of p-values from the null distribution
args = c(1:3)
args[1] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/MatrixEQTL/v8/all-by-all/'
args[2] = '/tigress/BEE/RNAseq/Output/trans-mapping/gtex/trans-eQTLs/MatrixEQTL/v8/all-by-all/'

summary_dir = args[1]
write_dir = args[2]

tissue_list = c('Adipose_Subcutaneous','Cells_EBV-transformed_lymphocytes','Skin_Sun_Exposed_Lower_leg','Testis','Thyroid','Whole_Blood')

library(dplyr)

trans_eqtl_list = list()

for (tissue in tissue_list) {
    print(tissue)
    list_files = list.files(path = paste0(summary_dir, tissue, '/summary/', sep=''))
    print(length(list_files))

    # Required data to aggregate
    cumulative_me_trans_total = vector('list', length(list_files))
    all_hist_vals_trans_total = rep(0,100)
    all_n_tests_trans_total = 0
    snps_tested_total = 0

    i = 1
    for (f in list_files) {
        load(paste0(summary_dir, tissue, '/summary/', f))
        
        cumulative_me_trans_total[[i]] = me_trans
        all_hist_vals_trans_total = all_hist_vals_trans_total + hist_counts
        all_n_tests_trans_total = all_n_tests_trans_total + total_ntests
        snps_tested_total = snps_tested_total + snps_tested
        i = i+1
    }

    cumulative_me_trans_total = do.call("rbind", cumulative_me_trans_total)
    print(dim(cumulative_me_trans_total))

    cumulative_me_trans_total = cumulative_me_trans_total[order(cumulative_me_trans_total$pvalue),]
    # BH procedure for SNP-gene FDR:
    cumulative_me_trans_total$FDR = p.adjust(cumulative_me_trans_total$pvalue, method = 'BH', n = all_n_tests_trans_total)

    n_eqtls_5 = sum(cumulative_me_trans_total$FDR <= 0.05)
    n_eqtls_10 = sum(cumulative_me_trans_total$FDR <= 0.10)

    print(n_eqtls_5)
    print(n_eqtls_10)

    print(length(unique(cumulative_me_trans_total$gene[1:n_eqtls_5])))
    print(length(unique(cumulative_me_trans_total$gene[1:n_eqtls_10])))
    print(length(unique(cumulative_me_trans_total$snps[1:n_eqtls_5])))
    print(length(unique(cumulative_me_trans_total$snps[1:n_eqtls_10])))

    trans_eqtl_list[[tissue]] = cumulative_me_trans_total[1:n_eqtls_5,]
    trans_eqtl_list[[tissue]]$tissue = tissue
}

trans_eqtls = do.call("rbind", trans_eqtl_list)

write.table(trans_eqtls, file = paste0(write_dir, 'v8_trans_eqtls_FDR_0.05.txt') sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)}
