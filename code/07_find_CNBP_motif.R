library(Biostrings)
library(dplyr)
scov2 <- as.character(paste0(read.table("../data/SARS-CoV-2.fa", skip = 1, header = FALSE)[,1], collapse = ""))
SCOV2 <- DNAString(scov2)
mp <- matchPattern("TGGAGNW", SCOV2, fixed = FALSE)
mp_RC <- matchPattern("TGGAGNW", reverseComplement(SCOV2), fixed = FALSE)
length(mp)
length(mp_RC)

perm <- sapply(1:1000, function(i){
  set.seed(i)
  SCOV2perm <- DNAString(paste0(sample(strsplit(scov2, "")[[1]]), collapse = ""))
  mp <- matchPattern("TGGAGNW", SCOV2perm, fixed = FALSE)
  mp_RC <- matchPattern("TGGAGNW", reverseComplement(SCOV2perm), fixed = FALSE)
  length(mp) 
})
z = (length(mp)-mean(perm))/sd(perm)

2*pnorm(-abs(z))
