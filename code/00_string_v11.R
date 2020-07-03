library(data.table)
library(dplyr)
"%ni%" <- Negate("%in%")

viral_genes <- c("N", "S", "3a", "M", "ORF9b", "PP1AB")

aliases <- fread("../string_v11/9606.protein.aliases.v11.0.txt.gz", skip = 1, header = FALSE, col.names = c("STRINGid", "Alias", "source"))[,c(1,2)]
string_v11 <- fread("../string_v11/9606.protein.links.v11.0.txt.gz")
gene2string <- (aliases$STRINGid); names(gene2string) <- (aliases$Alias)


pull_string_v11 <- function(genes_keep, min_score = 0){
  
  string_genes_keep <- gene2string[genes_keep]
  string2gene <- names(string_genes_keep); names(string2gene) <- unname(string_genes_keep)
  #print(paste0("WARNING: ", paste(genes_keep[is.na(string2gene)], collapse = ","), " not found"))

    # Now filter
  data <- string_v11 %>% filter(combined_score > min_score) %>%
    filter(protein1 %in% unname(string_genes_keep) & protein2 %in% unname(string_genes_keep))
  
  data$node1 <- unname(string2gene[as.character(data$protein1)])
  data$node2 <- unname(string2gene[as.character(data$protein2)])
  data[complete.cases(data),c("node1", "node2", "combined_score")]
  
}