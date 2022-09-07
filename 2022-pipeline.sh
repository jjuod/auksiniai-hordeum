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

# 4. map onto ref genome
#cd /mnt/quick/julius-temp
#wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/hordeum_vulgare_goldenpromise/dna_index/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa.gz
#bwa index -a bwtsw Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa.gz
cd ~/Documents/mieziai/2022/mapped/
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/AII_frw.fq.gz ../trim/AII_rev.fq.gz 2> bwa_log.txt | samtools sort -@2 -o AII_HVgp.bam
#bwa mem -t 9 /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	../trim/tw_frw.fq.gz ../trim/tw_rev.fq.gz 2> bwa_log_tw.txt | samtools sort -@2 -o tw_HVgp.bam
#samtools stats -in AII_HVgp.bam > map_stats_AII_HVgp.txt
#samtools stats -in tw_HVgp.bam > map_stats_tw_HVgp.txt
#plot-bamstats map_stats_AII_HVgp.txt -p plots/p_AII
#plot-bamstats map_stats_tw_HVgp.txt -p plots/p_tw
