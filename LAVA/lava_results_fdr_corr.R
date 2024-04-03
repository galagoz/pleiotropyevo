
lava_results = read.table("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/LAVA/D-RI_vs_wmTract_fromYasmina.txt",
                          header = T, sep = "\t")
##########################################################
# FDR
# remove duplicate rows (if there are any) before FDR
# because number of lines will be used as the number
# of independent tests by p.adjust
lava_results$fdr = p.adjust(lava_results$p, method = "fdr")
sig_lava_results = lava_results[lava_results$fdr<0.05,]

write.table(lava_results,
            "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/LAVA/lava_FDRall_results.csv", row.names = F, quote = F)