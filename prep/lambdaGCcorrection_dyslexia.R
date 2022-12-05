library(ggplot2)
library(data.table)
library(tidyverse)
library(cowplot)

# 1) Dyslexia

dys_ss = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2.dat")

chisq = qchisq(dys_ss$pvalue,1,lower.tail=FALSE)
newchisq = chisq/1.5733 # dylexia paper says lambdaGC is 1.18, but LDSC outputs 1.5733.
dys_ss$pvalue = pchisq(newchisq, df=1,lower.tail=FALSE)

fwrite(dys_ss, "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/dyslexia.filtered.2_lambdaGCcorrected.dat", quote = F, sep = "\t", row.names = F, na="NA")

dys_ss_hits <- filter(dys_ss,
                      dys_ss$rsid %in% paper_Pvals$SNP)

paper_Pvals$correctedPval = dys_ss_hits$pval_lambdaGC_corrected
############################
# raw p plot

p1 = ggplot(paper_Pvals, aes(x=Published_Pval, y=correctedPval)) + 
  geom_point()+
  geom_abline()+
  ylim(c(min(paper_Pvals$Published_Pval),max(paper_Pvals$correctedPval)))+
  xlim(c(min(paper_Pvals$Published_Pval),max(paper_Pvals$correctedPval)))+
  xlab("Published") + ylab("LambdaGC-corrected")

p2 = ggplot(paper_Pvals, aes(x=Published_Pval, y=Our_Pval)) + 
  geom_point()+
  geom_abline()+
  ylim(c(min(paper_Pvals$Published_Pval),max(paper_Pvals$Our_Pval)))+
  xlim(c(min(paper_Pvals$Published_Pval),max(paper_Pvals$Our_Pval)))+
  xlab("Published") + ylab("Uncorrected")

plot_grid(p2, p1)
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/figures/dyslexia_raw_pvals.pdf",
      width = 25, height = 12, units = "cm")
############################
# -log10(p) plot

p1_log = ggplot(paper_Pvals, aes(x=-log10(Published_Pval), y=-log10(correctedPval))) + 
  geom_point()+
  geom_abline()+
  ylim(c(min(-log10(paper_Pvals$correctedPval)),max(-log10(paper_Pvals$Published_Pval))))+
  xlim(c(min(-log10(paper_Pvals$correctedPval)),max(-log10(paper_Pvals$Published_Pval))))+
  xlab("Published") + ylab("LambdaGC-corrected")

p2_log = ggplot(paper_Pvals, aes(x=-log10(Published_Pval), y=-log10(Our_Pval))) + 
  geom_point()+
  geom_abline()+
  ylim(c(min(-log10(paper_Pvals$Published_Pval)),max(-log10(paper_Pvals$Our_Pval))))+
  xlim(c(min(-log10(paper_Pvals$Published_Pval)),max(-log10(paper_Pvals$Our_Pval))))+
  xlab("Published") + ylab("Uncorrected")

plot_grid(p2_log, p1_log)
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/figures/dyslexia_log_pvals.pdf",
       width = 25, height = 12, units = "cm")

# 2) Rhythm impairment

header = read.table("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythym_impairment_raw_sumstat.txt",
                     header = TRUE, nrow = 1)
rhyimp_ss = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythym_impairment_raw_sumstat.txt", 
                skip=1, header = FALSE)
setnames(rhyimp_ss, colnames(header))

chisq = qchisq(rhyimp_ss$pvalue,1,lower.tail=FALSE)
newchisq = chisq/1.4745
rhyimp_ss$pvalue = pchisq(newchisq, df=1,lower.tail=FALSE)

fwrite(rhyimp_ss, "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/data/rhythm_impairment_lambdaGCcorrected.dat", quote = F, sep = "\t", row.names = F, na="NA")

# 3) Genomic SEM CPM results

gensem_ss = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_dys_rhyimp.tab")

chisq = qchisq(gensem_ss$Pval_Estimate,1,lower.tail=FALSE)
newchisq = chisq/1.671
gensem_ss$Pval_Estimate = pchisq(newchisq, df=1,lower.tail=FALSE)

fwrite(gensem_ss, "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_dys_rhyimp_lambdaGCcorrected.tab", quote = F, sep = "\t", row.names = F, na="NA")

N_hat_F1<-mean(1/((2*gensem_ss$MAF*(1-gensem_ss$MAF))*gensem_ss$SE^2), na.rm = T)
gensem_ss$N = round(N_hat_F1)

fwrite(gensem_ss[,c(1,4,5,6,14,15,22)], file = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/GenomicSEM_multivarGWAS_dys_rhyimp_lambdaGCcorrected_forMunge.tab", 
       sep = "\t",  row.names = FALSE, col.names = TRUE)