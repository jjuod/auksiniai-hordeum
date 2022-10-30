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
# NOTE: HORVU.MOREX.r3.6HG0554760 not found
#grep -f genenames_melatonin.txt ../../refs/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.54.gff3 | awk '$3=="gene"{split($9,a,";"); print $1":"$4"-"$5, a[1]}' 
#samtools faidx -r <(awk '{print $1}' regions_melatonin.txt) ${REFM3} | awk '$0~/^>/{getline g < "regions_melatonin.txt"; print $0, g; next} {print}' > geneseqs_melatonin.txt

#blastn -query geneseqs_melatonin.txt -subject ${REFGP} -out blastres_regions_melatonin.txt -outfmt "6 qseqid qlen qstart qend sseqid sstart send evalue length pident"

