# ========================================================
# Make tables of LDSC partitioned heritability results

# The results files (.results) from the LDSC run should be organized
# by annotation, with separate directories for each annotation. 

# Updated for the nonGC version

# ========================================================

library(tidyverse)

options(stringsAsFactors=FALSE)

annots = dir(path = "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results/partherit",
                   full.names = F,
                   recursive = F, pattern = "results.results")

partheritresults = data.frame(Category = character(0),
                              Prop._SNPs = numeric(0),
                              Prop._h2 = numeric(0),
                              Prop._h2_std_error = numeric(0),
                              Enrichment = numeric(0),
                              Enrichment_std_error = numeric(0),
                              Enrichment_p = numeric(0),
                              Coefficient = numeric(0),
                              Coefficient_std_error = numeric(0),
                              Coefficient_z.score = numeric(0),
                              Annotation = character(0))

#i=2
for (i in 1:length(annots)){
      
      tmp_annot = strsplit(annots[i], split = "\\.")[[1]][2]
      tmp_annot
      
      file = Sys.glob(path = paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results/partherit/globularity.", 
                                     tmp_annot, ".results.results"))
      
      results = read.table(file,header=TRUE)
      results$annot = tmp_annot
      partheritresults[i,] = results[1,]
      
}

      partheritresults = partheritresults %>% 
        mutate(fdr = p.adjust(Enrichment_p, method = "fdr")) %>% # correcting for 3 tests
        ungroup()
      
      partheritresults$annot.p = if_else(partheritresults$fdr < 0.05, as.character(round(partheritresults$fdr, digits = 4)), "")
      
      partheritresults$significant = if_else(partheritresults$fdr < 0.05, "Yes", "")
      
      write.table(partheritresults, 
                  paste0("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results/partherit/globularity_partherit_results_FDR3.txt"),
                         sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)