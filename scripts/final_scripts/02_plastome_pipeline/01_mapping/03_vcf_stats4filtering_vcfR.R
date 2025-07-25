library(vcfR)
vcf <- read.vcfR("./data/Lb_plastomes_raw/Consensus/aln/all/all_plastomes.vcf.gz", verbose = FALSE )
dna <- ape::read.dna("./data/Lb_plastomes_raw/refs/NC_036356_1_Lus_chloroplast_ref.fasta", format = "fasta")
gff <- read.table("./data/Lb_plastomes_raw/refs/annotation.gff3", sep="\t", quote="")

chrom <- create.chromR(name = 'Chloroplast Genome', vcf=vcf, seq=dna, ann=gff)
chrom <- proc.chromR(chrom)
plot(chrom)
chromoqc(chrom)

## Based on plot
#10000 < DP < 30000
#MQ = 40
#QUAL

chrom <- masker(chrom, min_QUAL = 20, min_DP = 10000, max_DP = 40000, min_MQ = 30)
plot(chrom)
chromoqc(chrom)
