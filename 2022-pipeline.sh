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
cd ~/Documents/mieziai/2022/trim/
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -j 8 \
	-o tw1_frw.fq.gz \
	-p tw1_rev.fq.gz \
	../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_1.fq.gz \
	../raw/tw_EKDN220023658-1A_H33FGDSX5_L2_2.fq.gz &> errors_tw1.log
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -j 8 \
	-o tw2_frw.fq.gz \
	-p tw2_rev.fq.gz \
	../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_1.fq.gz \
	../raw/tw_EKDN220023658-1A_HYV73DSX3_L3_2.fq.gz &> errors_tw2.log
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -j 8 \
	-o tw3_frw.fq.gz \
	-p tw3_rev.fq.gz \
	../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_1.fq.gz \
	../raw/tw_EKDN220023658-1A_HYV73DSX3_L4_2.fq.gz &> errors_tw3.log
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
	-O 5 -q 20,20 -m 30 -j 8 \
	-o AII_frw.fq.gz \
	-p AII_rev.fq.gz \
	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_1.fq.gz \
	../raw/AII_EKDN220023659-1A_HTMGWDSX3_L1_2.fq.gz &> errors_AII.log
echo "Adapter trimming complete"
