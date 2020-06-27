library(data.table)
library(dplyr)

aliases <- fread("../string_v11/9606.protein.aliases.v11.0.txt.gz", skip = 1, header = FALSE, col.names = c("STRINGid", "Alias", "source"))[,c(1,2)]
data <- fread("../string_v11/9606.protein.links.v11.0.txt.gz")
trans_vec <- 
