# java -jar /mnt/hdd/soft/picard-2.19.0/picard.jar AddOrReplaceReadGroups \
#     I=ES_11_aligned_sorted.bam \
#     O=ES_11_aligned_sorted2.bam \
#     RGLB=ES RGPL=illumina RGPU=unit1 RGSM=11 \
#     VALIDATION_STRINGENCY=SILENT
# java -jar /mnt/hdd/soft/picard-2.19.0/picard.jar AddOrReplaceReadGroups \
#     I=ES_12_aligned_sorted.bam \
#     O=ES_12_aligned_sorted2.bam \
#     RGLB=ES RGPL=illumina RGPU=unit1 RGSM=12 \
#     VALIDATION_STRINGENCY=SILENT
# java -jar /mnt/hdd/soft/picard-2.19.0/picard.jar AddOrReplaceReadGroups \
#     I=PS_11_aligned_sorted.bam \
#     O=PS_11_aligned_sorted2.bam \
#     RGLB=PS RGPL=illumina RGPU=unit1 RGSM=11 \
#     VALIDATION_STRINGENCY=SILENT
# java -jar /mnt/hdd/soft/picard-2.19.0/picard.jar AddOrReplaceReadGroups \
#     I=PS_12_aligned_sorted.bam \
#     O=PS_12_aligned_sorted2.bam \
#     RGLB=PS RGPL=illumina RGPU=unit1 RGSM=12 \
#     VALIDATION_STRINGENCY=SILENT

# seq 1 7 | parallel -j 2 \
# 	/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java \
# 	-jar /mnt/hdd/soft/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar \
# 	HaplotypeCaller -R Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa \
# 	-I ES_11_aligned_sorted.bam \
# 	-I PS_11_aligned_sorted.bam \
# 	-L chr{}H \
# 	-O calls/hc_11_chr{}.vcf

# seq 1 5 | parallel -j 5 \
# 	/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java \
# 	-jar /mnt/hdd/soft/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar \
# 	HaplotypeCaller -R Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa \
# 	-I ES_12_aligned_sorted.bam \
# 	-I PS_12_aligned_sorted.bam \
# 	-L chr{}H \
# 	-O calls/hc_12_chr{}.vcf

# replace chr names (dumb)
# for c in {1..7}
# do
#     sed 's/^chr\([0-7]\)H/\1H/' calls/hc_11_chr${c}.vcf > calls/renamed_11_chr${c}.vcf
#     sed 's/^chr\([0-7]\)H/\1H/' calls/hc_12_chr${c}.vcf > calls/renamed_12_chr${c}.vcf
# done

# Using snpEff
# seq 1 7 | parallel -j 7 'java -jar /mnt/hdd/soft/snpEff/snpEff.jar Hordeum_vulgare calls/renamed_11_chr{}.vcf -v -s calls/sesum_11_chr{}.html > calls/snpeff_11_chr{}.vcf 2> calls/seerr_11_chr{}.txt'
# seq 1 7 | parallel -j 7 'java -jar /mnt/hdd/soft/snpEff/snpEff.jar Hordeum_vulgare calls/renamed_12_chr{}.vcf -v -s calls/sesum_12_chr{}.html > calls/snpeff_12_chr{}.vcf 2> calls/seerr_12_chr{}.txt'

# bgzip + tabix
# for c in {4..7}
# do
# 	bgzip calls/snpeff_11_chr${c}.vcf
# 	tabix -C calls/snpeff_11_chr${c}.vcf.gz
# 	bgzip calls/snpeff_12_chr${c}.vcf
# 	tabix -C calls/snpeff_12_chr${c}.vcf.gz
# done
#bcftools concat calls/snpeff_11_chr*.vcf.gz -Oz -o calls/snpeff_11_all.vcf.gz
#bcftools concat calls/snpeff_12_chr*.vcf.gz -Oz -o calls/snpeff_12_all.vcf.gz

## CALLING VIA MPILEUP
# seq 1 7 | parallel -j 4 'samtools mpileup \
#     -r chr{}H -u -f ../Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa \
#     ../ES_11_aligned_sorted.bam ../ES_12_aligned_sorted.bam \
#     ../PS_11_aligned_sorted.bam ../PS_12_aligned_sorted.bam |\
#     bcftools call -m -v -Oz -o called_chr{}.vcf.gz'
# 

/usr/lib/jvm/java-1.8.0-openjdk-amd64/bin/java \
	-jar /mnt/hdd/soft/gatk-4.1.1.0/gatk-package-4.1.1.0-local.jar \
	VariantFiltration \
	-R Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa \
	-V pileup/snpeff_all.vcf.gz \
	--filter-expression "QD<2.0 || FS > 60.0 || MQ < 30.0" \
	--filter-name "QDFSMQ" \
	-O pileup/snpeff_gatk.vcf.gz
