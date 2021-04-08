computeMatrix scale-regions -R hg38_merged_5utr.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 100 --skipZeros --outFileName CNBP_5utr.tab -o ../output/CNBP_5utr.tsv.gz
computeMatrix scale-regions -R hg38_merged_5utr.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 100 --skipZeros --outFileName LARP1_5utr.tab -o ../output/LARP1_5utr.tsv.gz

computeMatrix scale-regions -R hg38_merged_cds.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 500 --skipZeros --outFileName CNBP_cds.tab -o ../output/CNBP_cds.tsv.gz
computeMatrix scale-regions -R hg38_merged_cds.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 500 --skipZeros --outFileName LARP1_cds.tab -o ../output/LARP1_cds.tsv.gz

computeMatrix scale-regions -R hg38_merged_3utr.bed -S CNBP_ns.bw -b 0 -a 0 --regionBodyLength 100 --skipZeros --outFileName CNBP_3utr.tab -o ../output/CNBP_3utr.tsv.gz
computeMatrix scale-regions -R hg38_merged_3utr.bed -S LARP1_ns.bw -b 0 -a 0 --regionBodyLength 100 --skipZeros --outFileName LARP1_3utr.tab -o ../output/LARP1_3utr.tsv.gz


