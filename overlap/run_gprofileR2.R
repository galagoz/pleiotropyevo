# Original script is taken from Jason Stein's ENIGMA-evo repository.
# Gokberk Alagoz - 01/10/2021

library(gprofiler2)
options(stringsAsFactors = FALSE)

#Load in global gene list

GlobalSAHGE = read.csv("/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results/VEP/glob_gene_list.txt")

GO <- gost(GlobalSAHGE$x,
    organism = "hsapiens",
    exclude_iea = TRUE,
    correction_method = "fdr",
    significant = TRUE,
    sources = c("GO", "KEGG", "REAC"),
    ordered_query = FALSE,
    numeric_ns = "ENTREZGENE_ACC",
    user_threshold = 0.05,
    domain_scope = "annotated"
  )

GO_df = as.data.frame(GO$result)
GO_df_2 <- apply(GO_df, 2, as.character) # "flatten" your list into a character vector, so that you can write it.

write.table(GO_df_2, 
            "/data/workspaces/lag/workspaces/lg-ukbiobank/projects/globularity/evolution/results/VEP/GO_results.txt", 
            row.names = F, quote = F)
