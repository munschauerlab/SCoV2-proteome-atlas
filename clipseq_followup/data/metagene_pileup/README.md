# How to reproduce metagene pileup matricies

1) Download GTF from Ensembel or wherever (we used hg38)
2) Parse GTF into features using this: https://github.com/stephenfloor/extract-transcript-regions
3) Use bedtools merge to get consolidated features. 
4) Use computeMatrix (deeptools) to get meta gene profiles-- see `l2.sh`
5) rearrange the files a bit to then run the `04_metagene_pileup.R`

