set -e

# Post-calling things:
# extracting target genes etc

REFGP="/mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna_sm.toplevel.fa"
REFM3="/mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa"

# (run analysis.Rmd first and upload regions_*.txt)
cd ~/Documents/mieziai/2022/genes/
#samtools faidx -r regions_mod.txt ${REFGP} > geneseqs_mod.txt
#samtools faidx -r regions_hi.txt ${REFGP} > geneseqs_hi.txt

#blastn -query geneseqs_hi.txt -subject ${REFM3} -out blastres_regions_hi.txt -outfmt "6 qseqid qlen qstart qend sseqid sstart send evalue length pident"
#blastn -query geneseqs_mod.txt -subject ${REFM3} -out blastres_regions_mod.txt -outfmt "6 qseqid qlen qstart qend sseqid sstart send evalue length pident"

# extract full gene region sequences
#samtools faidx -r regions_buf1kb_mod.txt ${REFGP} > geneseqs_buf1kb_mod.fa
#samtools faidx -r regions_buf1kb_hi.txt ${REFGP} > geneseqs_buf1kb_hi.fa
#samtools faidx geneseqs_buf1kb_hi.fa
#samtools faidx geneseqs_buf1kb_mod.fa

# Melatonin genes
# (extract their variants from each individual)
# NOTE: HORVU.MOREX.r3.6HG0554760 not found
cd ~/Documents/mieziai/2022/genes/
grep -f genenames_melatonin.txt ../../refs/Hordeum_vulgare_goldenpromise.GPv1.54.gff3 | awk -v OFS="\t" '$3=="gene"{split($9,a,";"); print $1, $4-1000, $5+1000, a[1]}' > regions_melatonin.txt

cd ~/Documents/mieziai/old/called/vcfs/
bcftools concat called_AII-11_chr*.vcf.gz -n -o called_AII-11.vcf.gz
bcftools concat called_tw-12_chr*.vcf.gz -n -o called_tw-12.vcf.gz
bcftools index called_AII-11.vcf.gz
bcftools index called_tw-12.vcf.gz
bcftools merge called_AII-11.vcf.gz called_tw-12.vcf.gz -0 -Oz -o called_both1112.vcf.gz 
bcftools index called_both1112.vcf.gz

bcftools view -R ~/Documents/mieziai/2022/genes/regions_melatonin.txt called_both1112.vcf.gz -Ou |\
	bcftools query -H -f '%CHROM %POS %REF %ALT %MQ0F %DP %DP4[ %GT]\n' > \
	called_melatonin_GP_both1112.csv

