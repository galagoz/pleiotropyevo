#!/usr/bin/env Rscript
#
# This script will...
#
# Gokberk Alagoz
# created on: 12.2022

#####
library(ggplot2)
library(tidyverse)
library(data.table)
library(ggalt)
library(viridis)
library(cowplot)
library(ggpubr)
if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
options(stringsAsFactors=FALSE)

#--------------------------------------------

input = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_forEVO.txt")
outDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/"

colnames(input) = c("CHR", "STR", "END", "BETA", "SE", "P", "Q", "Q_P", "phastCons", "phyloP")

# Significantly hetero SNPs
input_QPranked = input[order(Q_P),]
input_sig = input[input$Q_P<5e-8,]

input_Qranked = input[order(Q),]
input_sig = input[input$P<5e-8,]
input_nomsig = input[input$P<0.05,]
nrow(input_nomsig)

input_nomsig$mlog10P = -log(input_nomsig$P)
input_nomsig$mlog10Q_P = -log(input_nomsig$Q_P)

# read in clumped SNPs list
clumped = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/plink/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_GCcorr_withN.clumped.leadSNP_genomicCoordinates.tab",
                header = T)
# remove NA rows
clumped = drop_na(clumped)
# add "chr" to the CHR column
clumped$CHR = sub("^", "chr", clumped$CHR)
colnames(clumped) = c("CHR", "END", "SNP")
# intersect input_nomsig and clumped
input_clumped = merge(input_nomsig, clumped, by = c("CHR", "END"))

nrow(input_clumped)
nrow(clumped)

# tag 18 genome-wide risk loci (this list is taken from FUMA)
input_nomsig$riskloci = "no"
input_nomsig[input_nomsig$CHR=="chr3"&input_nomsig$END=="17394615",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr3"&input_nomsig$END=="81811784",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr3"&input_nomsig$END=="135773661",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr4"&input_nomsig$END=="92035550",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr4"&input_nomsig$END=="152286935",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr4"&input_nomsig$END=="176874330",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr5"&input_nomsig$END=="114286650",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr6"&input_nomsig$END=="108906200",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr7"&input_nomsig$END=="106863170",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr11"&input_nomsig$END=="30395895",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr11"&input_nomsig$END=="111916647",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr12"&input_nomsig$END=="60971306",]$riskloci = "yes" # this SNP does not have a phastCons score!
input_nomsig[input_nomsig$CHR=="chr13"&input_nomsig$END=="59565064",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr16"&input_nomsig$END=="30125840",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr17"&input_nomsig$END=="34908385",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr17"&input_nomsig$END=="44076063",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr20"&input_nomsig$END=="31093514",]$riskloci = "yes"
input_nomsig[input_nomsig$CHR=="chr20"&input_nomsig$END=="47616271",]$riskloci = "yes"

#--------------------------------------------
input_sig_hetero = input[input$Q_P<5e-8,]
nrow(input_sig_hetero) # 1142 sig hetero SNPs

input_sig_hetero_sig = input_sig[input_sig$Q_P<5e-8,]

#--------------------------------------------
input_hetero = input_sig[input_sig$Q>10,]

#--------------------------------------------
input_homo = input_sig[input_sig$Q<1,]

#--------------------------------------------
# plot
plot(input_hetero$Q, input_hetero$phyloP)
plot(input_homo$BETA, input_homo$phyloP)
median(input_homo$phyloP)
median(input_homo$phastCons)

plot(input_hetero$Q, input_hetero$phastCons)
plot(input_homo$Q, input_homo$phastCons)
plot(input_hetero$BETA, input_hetero$phyloP)
plot(input_hetero$BETA, input_hetero$phastCons)
median(input_hetero$phyloP)
median(input_hetero$phastCons)

plot(input_sig$BETA, input_sig$phyloP)

# scatter - version 1
# phyloP vs. -log10(Q_P)
ggplot(NULL, aes(Q, phyloP)) +
       geom_point(data = input_sig) +
       geom_point(data = input_homo, color = "green") +
       geom_point(data = input_hetero, color = "red")

ggplot(NULL, aes(Q, phyloP)) +
  geom_point(data = input_sig, aes(colour = P))

p1 = ggplot(input_nomsig, aes(mlog10Q_P, phyloP)) +
  geom_point(aes(colour = mlog10P)) +
  geom_point(data = input_nomsig[input_nomsig$riskloci=="yes",], color = "red") +
  #geom_point(data = input_nomsig[input_nomsig$P<5e-9,], color = "red") +
  geom_vline(aes(xintercept = -log10(5e-08)), size = 1, lty = 2, color = "red") +
  border()

#ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/Qpvalues_vs_phyloP_CPM_FUMAriskLociMarked.png",
#       p1, width = 5, height = 5, dpi = 600)

# chr13:59565064 SNP has the lowest phyloP (strongest accelated selection)
# look further into this SNP.

# scatter - version 2
# phyloP vs. -log10(GWAS_P)

p2 = ggplot(input_nomsig, aes(mlog10P, phyloP)) +
  geom_point(aes(colour = mlog10Q_P)) +
  geom_point(data = input_nomsig[input_nomsig$riskloci=="yes",], color = "red") +
  geom_vline(aes(xintercept = -log10(5e-08)), size = 1, lty = 2, color = "red") +
  border()

#ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/GWASpvalues_vs_phyloP_CPM_FUMAriskLociMarked_v2.png",
#       p2, width = 5, height = 5, dpi = 600)

# scatter 3
# phastCons vs. -log10(GWAS_P)

p3 = ggplot(input_nomsig, aes(mlog10P, phastCons)) +
  geom_point(aes(colour = mlog10Q_P)) +
  geom_point(data = input_nomsig[input_nomsig$riskloci=="yes",], color = "red") +
  geom_vline(aes(xintercept = -log10(5e-08)), size = 1, lty = 2, color = "red") +
  border()

ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/GWASpvalues_vs_phastCons_CPM_FUMAriskLociMarked.png",
       p3, width = 5, height = 5, dpi = 600)

# extremely conserved one is chr11:111916647

# scatter 4 - with clumped SNPs
# phastCons vs. -log10(GWAS_P)

p4 = ggplot(input_clumped, aes(mlog10P, phastCons)) +
  geom_point(aes(colour = mlog10Q_P)) +
  geom_point(data = input_nomsig[input_nomsig$riskloci=="yes",], color = "red") +
  geom_vline(aes(xintercept = -log10(5e-08)), size = 1, lty = 2, color = "red") +
  border()

ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/GWASpvalues_vs_phastCons_CPM_FUMAriskLociMarked_clumped.pdf",
       p4, width = 5, height = 5, dpi = 600)

# extremely conserved one is chr11:111916647

p5 = ggplot(input_clumped, aes(phastCons, mlog10P)) +
  geom_point(colour = "#214D9D") +
  geom_point(data = input_nomsig[input_nomsig$riskloci=="yes",], color = "red") +
  geom_hline(aes(yintercept = -log10(5e-08)), size = 1, lty = 2, color = "#AA1926") +
  border() + theme(legend.position = "non")

ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/GWASpvalues_vs_phastCons_CPM_FUMAriskLociMarked_clumped_v2_woLegend.pdf",
       p5, width = 5, height = 5, dpi = 600)

#--------------------------------------------
#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP.pdf", width = 5, height = 5)
#ggplot(input, aes(x=Q, y=phyloP)) + 
#  geom_point()+
#  xlab("Q") + ylab("phyloP")
#dev.off()

ggscatter(input, x = "Q", y = "phyloP", 
          add = "reg.line", add.params = list(color = "red"),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "right", label.y.npc = "bottom"),
          xlab = "Q-scores", ylab = "phyloP")
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_v2.png", width = 5, height = 5, dpi = 300)``

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons.pdf", width = 5, height = 5)
#ggplot(input, aes(x=Q, y=phastCons)) + 
#  geom_point()+
#  xlab("Q") + ylab("phastCons")
#dev.off()

ggscatter(input, x = "Q", y = "phastCons", 
          add = "reg.line", add.params = list(color = "red"),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "right", label.y.npc = "top"),
          xlab = "Q-scores", ylab = "phastCons")
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_v2.png", width = 5, height = 5, dpi = 300)

#--------------------------------------------

ggscatter(input_hetero, x = "Q", y = "phyloP", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = c(6e-08,-6),
          ylim = c(-6,2))

ggscatter(input_homo, x = "Q", y = "phyloP", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.coord = c(51,-4.8),
          ylim = c(-6,2))

p1 = ggplot(input_hetero, aes(x=Q, y=phyloP)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Significantly heterogeneous effect-SNPs", 
       subtitle = "n = 1142")+
  xlab("Q") + ylab("phyloP")+
  sm_statCorr(color = '#0f993d', corr_method = 'spearman',
              linetype = 'dashed')

p2 = ggplot(input_homo, aes(x=Q, y=phyloP)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1142")+
  xlab("Q") + ylab("phyloP")

p1_2 = ggplot(input_hetero, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Significantly heterogeneous effect-SNPs", 
       subtitle = "n = 1142")+
  xlab("Q") + ylab("phastCons")

p2_2 = ggplot(input_homo, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1142")+
  xlab("Q") + ylab("phastCons")

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2, p1, labels = c('A', 'B'), label_size = 12)
dev.off()

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2_2, p1_2, labels = c('A', 'B'), label_size = 12)
dev.off()
#--------------------------------------------

p1 = hist(input_hetero$Q)
p2 = hist(input_hetero$phyloP)
plot(p1, col=rgb(0,0,1,1/4))
plot(p2, col=rgb(1,0,0,1/4))
#--------------------------------------------
#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP.pdf", width = 5, height = 5)
#ggplot(input, aes(x=Q, y=phyloP)) + 
#  geom_point()+
#  xlab("Q") + ylab("phyloP")
#dev.off()

ggscatter(input, x = "Q", y = "phyloP", 
          add = "reg.line", conf.int = TRUE, add.params = list(color = "red"),
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "right", label.y.npc = "bottom"),
          xlab = "Q-scores", ylab = "phyloP")
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_v2.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons.pdf", width = 5, height = 5)
#ggplot(input, aes(x=Q, y=phastCons)) + 
#  geom_point()+
#  xlab("Q") + ylab("phastCons")
#dev.off()

ggscatter(input, x = "Q", y = "phastCons", 
          add = "reg.line", conf.int = TRUE, add.params = list(color = "red"),
          cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.sep = "\n", label.x.npc = "right", label.y.npc = "top"),                                
          cor.coef.coord = c(6e-08,-6), xlab = "Q-scores", ylab = "phastCons")
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_v2.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter.pdf", width = 5, height = 5)
smoothScatter(input$Q, input$phyloP, transformation = function(x) x ^ 0.4,
              colramp = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                           "#FCFF00", "#FF9400", "#FF3100")))
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter.png", width = 5, height = 5, dpi = 300)


#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_hexagon.pdf", width = 5, height = 5)
p1 = ggplot(input, aes(Q, phyloP)) + 
  geom_hex(shape=16, size=0.10, show.legend = FALSE, color="red") +
  geom_smooth(method='lm', formula= y~x)+
  scale_fill_viridis()
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter_v2.png",
       p1, width = 5, height = 5, dpi = 300)
cor.test(input$Q, input$phyloP, method=c("pearson"))

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_hexagon.pdf", width = 5, height = 5)
p2 = ggplot(input, aes(Q, phastCons)) + 
  geom_hex(shape=16, size=0.25, show.legend = FALSE, color="red") +
  geom_smooth(method='lm', formula= y~x)+
  scale_fill_viridis()
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_smoothScatter_v2.png", 
       p2, width = 5, height = 5, dpi = 300)
cor.test(input$Q, input$phastCons, method=c("pearson"))


pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_phastCons.pdf", width = 10, height = 5)
plot_grid(p1, p2)
dev.off()

#--------------------------------------------
# make a LocusZoom plot of chr11:111916647 haplotype

input_chr11 = input[input$CHR=="chr11",]
input_rs10891314_0 = input_chr11[input_chr11$END>111600000,]
input_rs10891314 = input_rs10891314_0[input_rs10891314_0$END<112000000,]

locusZoom1 = ggplot(input_rs10891314, aes(END, -log10(P))) +
             geom_point(aes(colour = -log10(Q_P))) +
             border() + theme(legend.position = "none",
                              axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
locusZoom2 = ggplot(input_rs10891314, aes(END, phastCons)) +
             geom_line(aes(END, phastCons)) +
             border() + theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank())
locusZoom3 = ggplot(input_rs10891314, aes(END, phyloP)) +
             geom_line(aes(END, phyloP)) +
             border() + theme(axis.title.x=element_blank())

locusZoom_final = plot_grid(locusZoom1, locusZoom2, locusZoom3,
                            nrow = 3, ncol = 1, rel_heights = c(1, 0.3, 0.3),
                            align = "v")
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/rs10891314_locusZoom_Q_phastCons_phyloP.pdf",
       locusZoom_final,
       width = 5, height = 7, dpi = 600)

# plot locusZoom1 with legend
locusZoom1_2 = ggplot(input_rs10891314, aes(END, -log10(P))) +
  geom_point(aes(colour = -log10(Q_P))) +
  border() + theme(axis.title.x=element_blank(),
                   axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())
ggsave("/data/workspaces/lag/workspaces/lg-genlang/Working/23andMe/Dyslexia2/Evolution/dys_rhy_pleiotropy/results/gSEM/rs10891314_locusZoom_Q.pdf",
       locusZoom1_2,
       width = 5, height = 5, dpi = 600)
#--------------------------------------------

sessionInfo()
