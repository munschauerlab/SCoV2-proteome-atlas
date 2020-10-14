bsub -q big computeMatrix scale-regions -R hg38_merged_5utr.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 1000 --skipZeros --outFileName CNBP_5utr.tab -o CNBP_5utr.tsv.gz
bsub -q big computeMatrix scale-regions -R hg38_merged_5utr.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 1000 --skipZeros --outFileName LARP1_5utr.tab -o LARP1_5utr.tsv.gz

bsub -q big computeMatrix scale-regions -R hg38_merged_cds.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 5000 --skipZeros --outFileName CNBP_cds.tab -o CNBP_cds.tsv.gz
bsub -q big computeMatrix scale-regions -R hg38_merged_cds.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 5000 --skipZeros --outFileName LARP1_cds.tab -o LARP1_cds.tsv.gz

bsub -q big computeMatrix scale-regions -R hg38_merged_3utr.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 1000 --skipZeros --outFileName CNBP_3utr.tab -o CNBP_3utr.tsv.gz
bsub -q big computeMatrix scale-regions -R hg38_merged_3utr.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 1000 --skipZeros --outFileName LARP1_3utr.tab -o LARP1_3utr.tsv.gz


