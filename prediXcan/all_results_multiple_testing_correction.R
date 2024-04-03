library(dplyr)

twas_results = read.csv("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_all_results.csv", header = T) 
twas_results$pvalue = as.numeric(twas_results$pvalue)
twas_results = na.omit(twas_results)

nrow(twas_results)
# 193,541 gene-tissue pairs to be removed
length(unique(twas_results$gene))
# 20,167 unique genes

write.table(unique(twas_results$gene),
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_background.txt", row.names = F, quote = F)

# filter chr17 inversion region
chr17_inversion = read.table("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/chr17_inv_gns.txt", header = F)
chr17_inversion = as.vector(chr17_inversion$V1)

tobefiltered = twas_results %>% filter(gene %in% chr17_inversion)
nrow(tobefiltered)
# 636 gene-tissue pairs to be removed
length(unique(tobefiltered$gene))
# 62 genes to be removed

twas_results_filtered = twas_results %>% filter(!gene %in% tobefiltered$gene)
nrow(twas_results_filtered)
# 192,905 gene-tissue pairs left
length(unique(twas_results_filtered$gene_name))
# 20,102 unique genes

##########################################################
# Bonferroni
# Bonferroni Threshold = 0.05 / number of unique gene tested (20167) / number of tissues (14)
sig_twas_results = twas_results[twas_results$pvalue < 0.05/20167/14,]

nrow(sig_twas_results)
# 223 significant genes

write.table(sig_twas_results,
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_sig_results.csv", row.names = F, quote = F)

# 36 unique genes
length(unique(sig_twas_results$gene))
write.table(unique(sig_twas_results$gene),
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_unique_sig_genes.txt", row.names = F, quote = F)

##########################################################
# FDR
# remove duplicate rows (if there are any) before FDR
# because number of lines will be used as the number
# of independent tests by p.adjust
twas_results_filtered = twas_results_filtered[!duplicated(twas_results_filtered),]
twas_results_filtered$fdr = p.adjust(twas_results_filtered$pvalue, method = "fdr")
sig_twas_results = twas_results_filtered[twas_results_filtered$fdr<0.05,]

nrow(sig_twas_results)
# 1275 significant gene-tissue pairs
length(unique(sig_twas_results$gene))
# 315 unique genes

sig_twas_results = sig_twas_results[order(sig_twas_results$fdr),]

write.table(twas_results_filtered,
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_FDRall_results.csv", row.names = F, quote = F)
write.table(sig_twas_results,
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_FDRsig_results.csv", row.names = F, quote = F)
write.table(unique(sig_twas_results$gene),
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/TWAS/twas_unique_FDRsig_genes.txt", row.names = F, quote = F)
