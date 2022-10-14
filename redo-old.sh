set -e

# 1. QC: skipping, was produced last time

# 2. cutadapt trimming. already done with:
#cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#    -o ES_11_frw.fastq.gz \
#    -p ES_11_rev.fastq.gz \
#    -m 30 \
#    -j 3 ${RAW_FILES}
# 11: AII, 12: tw

# 3. Repeat QC after trimming - done before, AND merge batches
cd ~/Documents/mieziai/old/trim/
#zcat ES_11_frw.fastq.gz PS_11_frw.fastq.gz | gzip -c > AII-11_frw.fq.gz
#zcat ES_11_rev.fastq.gz PS_11_rev.fastq.gz | gzip -c > AII-11_rev.fq.gz
#zcat ES_12_frw.fastq.gz PS_12_frw.fastq.gz | gzip -c > tw-12_frw.fq.gz
#zcat ES_12_rev.fastq.gz PS_12_rev.fastq.gz | gzip -c > tw-12_rev.fq.gz
echo "post-trim QC complete"

# 4. get and index both references
# skipping, using same GPv1

# 5. map onto ref genome
cd ~/Documents/mieziai/old/mapped/
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/AII-11_frw.fq.gz ../trim/AII-11_rev.fq.gz \
#	2> bwa_log.txt | samtools fixmate -@2 - - | \
#	samtools sort -@4 -m 3500M -o /mnt/quick/julius-temp/mapping/AII-11_HVgp.bam
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/tw-12_frw.fq.gz ../trim/tw-12_rev.fq.gz \
#	2> bwa_log_tw.txt | samtools fixmate -@2 - - | \
#       samtools sort -@4 -m 3500M -o /mnt/quick/julius-temp/mapping/tw-12_HVgp.bam
#samtools stats -in /mnt/quick/julius-temp/mapping/AII-11_HVgp.bam > map-stats-AII-11-HVgp.txt
#samtools stats -in /mnt/quick/julius-temp/mapping/tw-12_HVgp.bam > map-stats-tw-12-HVgp.txt
#plot-bamstats map-stats-AII-11-HVgp.txt -p plots/old-AII
#plot-bamstats map-stats-tw-12-HVgp.txt -p plots/old-tw

echo "mapping complete"

# 6. remove duplicates, while doing appropriate sorts
# NOTE: NOT DONE b/c i forgot the -m flag on fixmate. Shouldn't be a problem as the number of duplicates was small.

# 7. Calling via mileup
cd /mnt/quick/julius-temp/mapping/

#samtools index -c -@8 AII-11_HVgp.bam
#samtools index -c -@8 tw-12_HVgp.bam

# call separately then
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r contig{} -u -f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#    AII-11_HVgp.bam |\
#    bcftools call -m -v -Oz -o ../called_AII-11_chr{}.vcf.gz'
#
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r contig{} -u -f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#    tw-12_HVgp.bam |\
#    bcftools call -m -v -Oz -o ../called_tw-12_chr{}.vcf.gz'

echo "mpileup complete"

## convert calls to text format
cd ~/Documents/mieziai/old/called/
#for i in {1..7}
#do
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_AII-11_HVgp_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called_AII-11_chr${i}.vcf.gz >> calls_AII-11_HVgp_chr${i}.txt
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_tw-12_HVgp_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called_tw-12_chr${i}.vcf.gz >> calls_tw-12_HVgp_chr${i}.txt
#
#	awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' calls_AII-11_HVgp_chr${i}.txt > callsf_AII-11_HVgp_chr${i}.txt
#	awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' calls_tw-12_HVgp_chr${i}.txt > callsf_tw-12_HVgp_chr${i}.txt
#done


# TESTING
# generate consensus:
# normalize vcf
cd ~/Documents/mieziai/old/called/
/mnt/hdd/soft/bcftools-1.16/bcftools norm /mnt/quick/julius-temp/called_AII-11_chr4.vcf.gz \
       	-f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
	-Oz -o normd_AII-11_chr4.vcf.gz
bcftools index normd_AII-11_chr4.vcf.gz

# extract 1chr reference
grep -n contig /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa
awk 'NR>=19795980 && NR<29362982' /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa > ../../refs/Hordeum_vulgare_goldenpromise.dna.chr4.fa
bwa index -a bwtsw ../../refs/Hordeum_vulgare_goldenpromise.dna.chr4.fa

# create consensus
/mnt/hdd/soft/bcftools-1.16/bcftools consensus -HA -f ../../refs/Hordeum_vulgare_goldenpromise.dna.chr4.fa normd_AII-11_chr4.vcf.gz -i'DP>5' > consensus_AII-11_chr4.fa 2> cons_AII-11_chr4.log
# Around 2k/1.1M are skipped, mostly indels overlapping other variants

cd ~/Documents/mieziai/2022/mapped/

# extract mapped fastqs
samtools view -@4 -f PROPER_PAIR -b AII_HVgp_nodups.bam contig4 \
	-o /mnt/quick/julius-temp/mapping/AII_HVgp_nodups_chr4.bam
samtools collate -@4 -f -u -O /mnt/quick/julius-temp/mapping/AII_HVgp_nodups_chr4.bam |\
       	samtools fastq -@4 -1 AII_HVgp_nodups_chr4_1.fastq -2 AII_HVgp_nodups_chr4_2.fastq
gzip AII_HVgp_nodups_chr4_1.fastq
gzip AII_HVgp_nodups_chr4_2.fastq

bwa mem -t 9 ~/Documents/mieziai/old/called/consensus_AII-11_chr4.fa \
	AII_HVgp_nodups_chr4_1.fastq.gz AII_HVgp_nodups_chr4_2.fastq.gz \
	2> bwa_cons_chr4_log.txt | samtools fixmate -@2 -m - - | \
	samtools sort -@4 -m 3500M -o /mnt/quick/julius-temp/mapping/AII_HVgp_cons_chr4.bam

