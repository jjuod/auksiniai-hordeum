set -e

# 1. QC check
# cd ~/Documents/mieziai/2022/qc/
# zcat ../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_1.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:AII_L1_1_qc
# zcat ../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_2.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:AII_L1_2_qc
# zcat ../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_1.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L2_1_qc
# zcat ../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_2.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L2_2_qc
# zcat ../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_1.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L3_1_qc
# zcat ../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_2.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L3_2_qc
# zcat ../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_1.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L4_1_qc
# zcat ../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_2.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_L4_2_qc
# echo "QC complete"

# 2. cutadapt trimming
# (adapters taken from the vendor qc report)
#cd ~/Documents/mieziai/2022/trim/
#cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
#	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
#	-O 5 -q 20,20 -m 30 -j 8 \
#	-o tw1_frw.fq.gz \
#	-p tw1_rev.fq.gz \
#	../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_1.fq.gz \
#	../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_2.fq.gz &> errors_tw1.log
#cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
#	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
#	-O 5 -q 20,20 -m 30 -j 8 \
#	-o tw2_frw.fq.gz \
#	-p tw2_rev.fq.gz \
#	../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_1.fq.gz \
#	../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_2.fq.gz &> errors_tw2.log
#cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
#	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
#	-O 5 -q 20,20 -m 30 -j 8 \
#	-o tw3_frw.fq.gz \
#	-p tw3_rev.fq.gz \
#	../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_1.fq.gz \
#	../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_2.fq.gz &> errors_tw3.log
#cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
#	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
#	-O 5 -q 20,20 -m 30 -j 8 \
#	-o AII_frw.fq.gz \
#	-p AII_rev.fq.gz \
#	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_1.fq.gz \
#	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_2.fq.gz &> errors_AII.log
#echo "Adapter trimming complete"

# 3. Repeat QC after trimming AND merge tw batches
# (this can also be done w/ cat directly, but for downstream safety:)
#zcat tw1_frw.fq.gz tw2_frw.fq.gz tw3_frw.fq.gz | gzip -c > tw_frw.fq.gz
#zcat tw1_rev.fq.gz tw2_rev.fq.gz tw3_rev.fq.gz | gzip -c > tw_rev.fq.gz
#
#zcat AII_frw.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:AII_frw_qc
#zcat AII_rev.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:AII_rev_qc
#zcat tw_frw.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_frw_qc
#zcat tw_rev.fq.gz | /mnt/hdd/soft/FastQC/fastqc stdin:tw_rev_qc
#echo "post-trim QC complete"

# 4. get and index both references
#cd /mnt/quick/julius-temp
#wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/hordeum_vulgare_goldenpromise/dna_index/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa.gz
# gunzip Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa.gz
#bwa index -a bwtsw Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa.gz

#wget https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/plants/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna_sm.toplevel.fa.gz
#mv Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna_sm.toplevel.fa.gz Hordeum_vulgare.MorexV3.dna_sm.fa.gz
#gunzip Hordeum_vulgare.MorexV3.dna_sm.fa.gz
#bwa index -a bwtsw Hordeum_vulgare.MorexV3.dna_sm.fa

# 5. map onto ref genome
#cd ~/Documents/mieziai/2022/mapped/
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/AII_frw.fq.gz ../trim/AII_rev.fq.gz 2> bwa_log.txt | samtools sort -@2 -o AII_HVgp.bam
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/tw_frw.fq.gz ../trim/tw_rev.fq.gz 2> bwa_log_tw.txt | samtools sort -@2 -o tw_HVgp.bam
#samtools stats -in AII_HVgp.bam > map_stats_AII_HVgp.txt
#samtools stats -in tw_HVgp.bam > map_stats_tw_HVgp.txt
#plot-bamstats map_stats_AII_HVgp.txt -p plots/p_AII
#plot-bamstats map_stats_tw_HVgp.txt -p plots/p_tw

#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#	../trim/AII_frw.fq.gz ../trim/AII_rev.fq.gz 2> bwa_log.txt | samtools fixmate -m -@2 - AII_M3.bam
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#	../trim/tw_frw.fq.gz ../trim/tw_rev.fq.gz 2> bwa_log_tw.txt | samtools fixmate -m -@2 - tw_M3.bam
#samtools stats -in AII_M3.bam > map-stats-AII-M3.txt
#samtools stats -in tw_M3.bam > map-stats-tw-M3.txt
#plot-bamstats map-stats-AII-M3.txt -p plots/M3-AII
#plot-bamstats map-stats-tw-M3.txt -p plots/M3-tw

#samtools sort -@8 -m 5000M -o /mnt/quick/julius-temp/AII_M3_sorted.bam AII_M3.bam
#mv /mnt/quick/julius-temp/AII_M3_sorted.bam AII_M3.bam
#samtools sort -@8 -m 5000M -o /mnt/quick/julius-temp/tw_M3_sorted.bam tw_M3.bam
#mv /mnt/quick/julius-temp/tw_M3_sorted.bam tw_M3.bam

# M3 seems to have similar amount of total mapped reads, but a lot more of them at MQ0.


# 6. remove duplicates, while doing appropriate sorts
#samtools sort -n -@8 -m 3000M AII_HVgp.bam | samtools fixmate -@3 -m - - | samtools sort -@8 -m 3000M | samtools markdup -@3 -r -s - AII_HVgp_nodups.bam &> dup_stats_AII_HVgp.txt
#samtools sort -n -@8 -m 3000M tw_HVgp.bam | samtools fixmate -@3 -m - - | samtools sort -@8 -m 3000M | samtools markdup -@3 -r -s - tw_HVgp_nodups.bam &> dup_stats_tw_HV_gp.txt

#samtools index -c -@6 AII_HVgp_nodups.bam
#samtools index -c -@6 tw_HVgp_nodups.bam
#samtools index -c -@6 AII_M3.bam
#samtools index -c -@6 tw_M3.bam

# 7. Calling via mileup

# call jointly? but this doesn't output depth per each sample, only GTs and joint depth!
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r contig{} -u -f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#    AII_HVgp_nodups.bam tw_HVgp_nodups.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/called22_chr{}.vcf.gz'

# call separately then
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r contig{} -u -f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#    AII_HVgp_nodups.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/called22_AII_chr{}.vcf.gz'
#
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r contig{} -u -f /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#    tw_HVgp_nodups.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/called22_tw_chr{}.vcf.gz'
#
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r {}H -u -f /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#    AII_M3.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/called22_M3_AII_chr{}.vcf.gz'
#seq 1 7 | parallel -j 7 'samtools mpileup \
#    -r {}H -u -f /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#    tw_M3.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/called22_M3_tw_chr{}.vcf.gz'

## convert calls to text format
cd ~/Documents/mieziai/2022/called/
#for i in {1..7}
#do
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_AII_HVgp_nodups_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called22_AII_chr${i}.vcf.gz >> calls_AII_HVgp_nodups_chr${i}.txt
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_tw_HVgp_nodups_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called22_tw_chr${i}.vcf.gz >> calls_tw_HVgp_nodups_chr${i}.txt
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_AII_M3_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called22_M3_AII_chr${i}.vcf.gz >> calls_AII_M3_chr${i}.txt
#	echo 'POS REF ALT DP MQ MQ0F DP4 GT' > calls_tw_M3_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		/mnt/quick/julius-temp/called22_M3_tw_chr${i}.vcf.gz >> calls_tw_M3_chr${i}.txt
#done

#for CHR in {1..7}
#do
#	awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' calls_AII_HVgp_nodups_chr${CHR}.txt > callsf_AII_HVgp_nodups_chr${CHR}.txt
#	awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' calls_tw_HVgp_nodups_chr${CHR}.txt > callsf_tw_HVgp_nodups_chr${CHR}.txt
#done

# 8. ANNOTATE VIA SNPEFF

cd ~/Documents/mieziai/2022/called/
#bcftools concat /mnt/quick/julius-temp/called22_tw_chr1.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr2.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr3.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr4.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr5.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr6.vcf.gz \
#	/mnt/quick/julius-temp/called22_tw_chr7.vcf.gz \
#	-Ou -n | bcftools view  --min-af 1 \
#	-Oz -o called22_tw_mono.vcf.gz
#java -jar /mnt/hdd/soft/snpEff/snpEff.jar Hordeum_vulgare_goldenpromise called22_tw_mono.vcf.gz | gzip > annot_tw_mono.vcf.gz
# apply some filtering and reforamt
#bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP\t%DP4\t%MQ0F\t%MQ\t%ANN\t[%GT]\n" annot_tw_mono.vcf.gz | awk -v FS="\t" 'BEGIN{print "CHR", "POS", "REF", "ALT", "DP", "DP4", "MAF", "MQ0F", "MQ", "ANN"} $5>10 && $8>20 && $10=="1\/1"{split($6, d, ","); ref=d[1]+d[2]; alt=d[3]+d[4];
#print substr($1, 7,7), $2, $3, $4, $5, ref+alt, alt/(ref+alt), $7, $8, $9}' > annot_tw_mono.txt



# TESTING:
# calc raw read depth
cd ~/Documents/mieziai/2022/depth/
#samtools depth -a ../mapped/AII_HVgp_nodups.bam -r contig4 | cut -f 3 | gzip -c > dp_AII_HVgp_nodups_chr4.txt.gz
#samtools depth -a ../mapped/tw_HVgp_nodups.bam -r contig4 | cut -f 3 | gzip -c > dp_tw_HVgp_nodups_chr4.txt.gz
#samtools depth -a ../mapped/AII_HVgp_nodups.bam -r contig3 | cut -f 3 | gzip -c > dp_AII_HVgp_nodups_chr3.txt.gz
#samtools depth -a ../mapped/tw_HVgp_nodups.bam -r contig3 | cut -f 3 | gzip -c > dp_tw_HVgp_nodups_chr3.txt.gz

# average over 1kb
#zcat dp_AII_HVgp_nodups_chr3.txt.gz | awk '{a=a+$1} NR%1000==0{print int(NR/1000), a; a=0}' > dpkb_AII_HVgp_nodups_chr3.txt
#zcat dp_tw_HVgp_nodups_chr3.txt.gz | awk '{a=a+$1} NR%1000==0{print int(NR/1000), a; a=0}' > dpkb_tw_HVgp_nodups_chr3.txt
#zcat dp_AII_HVgp_nodups_chr4.txt.gz | awk '{a=a+$1} NR%1000==0{print int(NR/1000), a; a=0}' > dpkb_AII_HVgp_nodups_chr4.txt
#zcat dp_tw_HVgp_nodups_chr4.txt.gz | awk '{a=a+$1} NR%1000==0{print int(NR/1000), a; a=0}' > dpkb_tw_HVgp_nodups_chr4.txt

# ALTERNATIVE REFERENCE
#cd ~/Documents/mieziai/2022/mapped/
#samtools stats -in /mnt/quick/julius-temp/Barke/twL4_barke.bam > map-stats-twL4-barke.txt
#plot-bamstats map-stats-twL4-barke.txt -p plots/barke-twL4
#samtools mpileup \
#    -u -f /mnt/quick/julius-temp/Barke/assembly2_WGSBarke_renamed_blastable_carma.fasta \
#    /mnt/quick/julius-temp/Barke/twL4_barke.bam |\
#    bcftools call -m -v -Oz -o /mnt/quick/julius-temp/Barke/called22_barke_twL4.vcf.gz
#cd ~/Documents/mieziai/2022/called/
#echo 'CHR POS REF ALT DP MQ MQ0F DP4 GT' > calls_twL4_barke.txt
#bcftools query -f '%CHROM %POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#	/mnt/quick/julius-temp/Barke/called22_barke_twL4.vcf.gz >> calls_twL4_barke.txt
#awk 'NR==1{print $1,$2,$3,$4,$5,$6,$7,"DP4","MAF",$9; next} {split($8,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, $7, rr+aa, aa/(rr+aa), $9}' calls_twL4_barke.txt > callsf_twL4_barke.txt

# NOT DONE:
#/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java \
#	-jar /mnt/hdd/soft/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar \
#	VariantFiltration \
#	-R Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa \
#	-V pileup/snpeff_all.vcf.gz \
#	--filter-expression "QD<2.0 || FS > 60.0 || MQ < 30.0" \
#	--filter-name "QDFSMQ" \
#	-O pileup/snpeff_gatk.vcf.gz

