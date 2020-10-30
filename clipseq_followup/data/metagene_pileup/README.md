# How to reproduce metagene pileup matricies

1) Download GTF from Ensembel or wherever (we used hg38)
2a) Parse GTF into features using this: https://github.com/stephenfloor/extract-transcript-regions
2b) Use bedtools merge to get consolidated features. 
3) Use computeMatrix (deeptools)-- see `l2.sh`
4) rearrange the fiels a bit to then run the `04_metagene_pileup.R`

