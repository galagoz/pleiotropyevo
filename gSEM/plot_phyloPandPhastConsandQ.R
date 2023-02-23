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
options(stringsAsFactors=FALSE)

#--------------------------------------------

input = fread("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/GenomicSEM_multivarGWAS_CPM_dys_rhyimp_forEVO.txt")
outDir = "/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/"

colnames(input) = c("CHR", "STR", "END", "BETA", "SE", "P", "Q", "Q_P", "phastCons", "phyloP")

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP.pdf", width = 5, height = 5)
ggplot(input, aes(x=Q, y=phyloP)) + 
  geom_point()+
  xlab("Q") + ylab("phyloP")
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons.pdf", width = 5, height = 5)
ggplot(input, aes(x=Q, y=phastCons)) + 
  geom_point()+
  xlab("Q") + ylab("phastCons")
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_hexagon.pdf", width = 5, height = 5)
ggplot(input, aes(Q, phyloP)) + 
  geom_hex(shape=16, size=0.10, show.legend = FALSE, color="red") +
  scale_fill_viridis()
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter.pdf", width = 5, height = 5)
smoothScatter(input$Q, input$phyloP, transformation = function(x) x ^ 0.4,
              colramp = colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                           "#FCFF00", "#FF9400", "#FF3100")))
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_smoothScatter.png", width = 5, height = 5, dpi = 300)

#pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_hexagon.pdf", width = 5, height = 5)
ggplot(input, aes(Q, phastCons)) + 
  geom_hex(shape=16, size=0.25, show.legend = FALSE, color="red") +
  scale_fill_viridis()
#dev.off()
ggsave("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_hexagon.png", width = 5, height = 5, dpi = 300)
#--------------------------------------------
input_hetero = input[input$Q_P<5e-8,]
nrow(input_hetero) # 1181 sig hetero SNPs

input_Qranked = input[order(Q),]
input_hetero=head(input_Qranked, n = 1181)
input_homo=tail(input_Qranked, n = 1181)

p1 = ggplot(input_hetero, aes(x=Q, y=phyloP)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Significantly heterogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phyloP")

p2 = ggplot(input_homo, aes(x=Q, y=phyloP)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phyloP")

p1_2 = ggplot(input_hetero, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Significantly heterogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phastCons")

p2_2 = ggplot(input_homo, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phastCons")

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2, p1, labels = c('A', 'B'), label_size = 12)
dev.off()

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2_2, p1_2, labels = c('A', 'B'), label_size = 12)
dev.off()
#--------------------------------------------

# FINISH THIS PART

p1 = ggscatter(input_hetero, x = "Q", y = "phyloP", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "pearson",
               xlab = "Q", ylab = "phyloP", title = "Significantly heterogeneous effect-SNPs")

p2 = ggplot(input_homo, aes(x=Q, y=phyloP)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phyloP")

p1_2 = ggplot(input_hetero, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Significantly heterogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phastCons")

p2_2 = ggplot(input_homo, aes(x=Q, y=phastCons)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)+
  labs(title = "Least homogeneous effect-SNPs", 
       subtitle = "n = 1181")+
  xlab("Q") + ylab("phastCons")

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phyloP_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2, p1, labels = c('A', 'B'), label_size = 12)
dev.off()

pdf("/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/results/gSEM/evo_results/Q_phastCons_heteroVShomo.pdf", width = 8, height = 5)
plot_grid(p2_2, p1_2, labels = c('A', 'B'), label_size = 12)
dev.off()
#--------------------------------------------

sessionInfo()
