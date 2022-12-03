set -e

M3REF=/mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa 

# 1. QC: skipping

# 2. cutadapt trimming: skipping

# 3. Repeat QC after trimming: skipping

# 4. get and index both references: skipping

# 5. map onto ref genome
cd ~/Documents/mieziai/old/mapped_M3/
#bwa mem -t 7 ${M3REF} \
#	../trim/AII-11_frw.fq.gz ../trim/AII-11_rev.fq.gz \
#	2> bwa_log_AII.txt | \
#	samtools sort -@4 -m 3500M -o AII-11_M3.bam
#bwa mem -t 9 ${M3REF} \
#	../trim/tw-12_frw.fq.gz ../trim/tw-12_rev.fq.gz \
#	2> bwa_log_tw.txt | \
#       samtools sort -@4 -m 3500M -o tw-12_M3.bam
#samtools stats -in AII-11_M3.bam > map-stats-AII-11-M3.txt
#samtools stats -in tw-12_M3.bam > map-stats-tw-12-M3.txt
#plot-bamstats map-stats-AII-11-M3.txt -p plots/old-AII
#plot-bamstats map-stats-tw-12-M3.txt -p plots/old-tw

#samtools index -c -@8 AII-11_M3.bam
#samtools index -c -@8 tw-12_M3.bam

echo "mapping complete"

# 6. remove duplicates, while doing appropriate sorts
# NOTE: not done, the duplicate number was small on these.

# 7. Calling via mileup
cd ~/Documents/mieziai/old/mapped_M3/

# call separately then
#seq 1 7 | parallel -j 7 'bcftools mpileup \
#    -r {}H -Ou -f /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#    AII-11_M3.bam |\
#    bcftools call -m -v -Oz -o ../called_M3/vcfs/called_M3_AII-11_chr{}.vcf.gz'
#
#seq 1 7 | parallel -j 7 'bcftools mpileup \
#    -r {}H -Ou -f /mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa \
#    tw-12_M3.bam |\
#    bcftools call -m -v -Oz -o ../called_M3/vcfs/called_M3_tw-12_chr{}.vcf.gz'

echo "mpileup complete"

## convert calls to text format
cd ~/Documents/mieziai/old/called_M3/
#for i in {1..7}
#do
#	echo 'POS REF ALT DP MQ MQ0F DP4 MAF GT' > callsf_AII-11_M3_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		vcfs/called_M3_AII-11_chr${i}.vcf.gz |\
#		awk '{split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' >> callsf_AII-11_M3_chr${i}.txt
#
#	echo 'POS REF ALT DP MQ MQ0F DP4 MAF GT' > callsf_tw-12_M3_chr${i}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		vcfs/called_M3_tw-12_chr${i}.vcf.gz |\
#		awk '{split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' >> callsf_tw-12_M3_chr${i}.txt
#done

# 8. GENERATING AII-11 CONSENSUS
cd ~/Documents/mieziai/old/called_M3/
echo "generating consensus from old samples" 
#for CHR in {1..7}
#do
#	# normalize vcf to clean up indels
#	/mnt/hdd/soft/bcftools-1.16/bcftools norm vcfs/called_M3_AII-11_chr${CHR}.vcf.gz \
#	       	-f ${M3REF} -Oz -o normd_M3_AII-11_chr${CHR}.vcf.gz
#	bcftools index normd_M3_AII-11_chr${CHR}.vcf.gz
#	
#	# extract each chr references
#	awk '$0~/^>'${CHR}'H/{p=1;print;next} p==1{if($0~/^>/){exit}; print}' ${M3REF} > ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa
#	bwa index -a bwtsw ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa
#	
#	# create consensus, taking alternate allele in all variable positions
#	/mnt/hdd/soft/bcftools-1.16/bcftools consensus -HA \
#		-f ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa \
#		normd_M3_AII-11_chr${CHR}.vcf.gz -i'DP>5' > consensus_M3_AII-11_chr${CHR}.fa 2> cons_M3_AII-11_chr${CHR}.log
#	# Around 2k/1.1M per chr are skipped, mostly indels overlapping other variants.
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa 
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa.sa
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa.pac 
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa.amb
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa.ann 
#	rm ../../refs/Hordeum_vulgare_MorexV3.dna.chr${CHR}.fa.bwt 
#	bwa index -a bwtsw consensus_M3_AII-11_chr${CHR}.fa
#done

# 9. MAP NEW POOLS FROM THE CONSENSUS

cd ~/Documents/mieziai/2022/mapped/
echo "mapping the separate pools"

#for CHR in {1..7}
#do
#	# extract (properly) mapped fastqs
#	samtools view -@4 -f PROPER_PAIR -b AII_M3.bam ${CHR}H \
#	        -o /mnt/quick/julius-temp/mapping/AII_M3_chr${CHR}.bam
#	samtools collate -@4 -f -u -O /mnt/quick/julius-temp/mapping/AII_M3_chr${CHR}.bam |\
#	       	samtools fastq -@4 -1 AII_M3_chr${CHR}_1.fastq \
#		-2 AII_M3_chr${CHR}_2.fastq
#	gzip AII_M3_chr${CHR}_1.fastq
#	gzip AII_M3_chr${CHR}_2.fastq
#	rm /mnt/quick/julius-temp/mapping/AII_M3_chr${CHR}.bam
#
#	samtools view -@4 -f PROPER_PAIR -b tw_M3.bam ${CHR}H \
#		-o /mnt/quick/julius-temp/mapping/tw_M3_chr${CHR}.bam
#	samtools collate -@4 -f -u -O /mnt/quick/julius-temp/mapping/tw_M3_chr${CHR}.bam |\
#	       	samtools fastq -@4 -1 tw_M3_chr${CHR}_1.fastq \
#		-2 tw_M3_chr${CHR}_2.fastq
#	gzip tw_M3_chr${CHR}_1.fastq
#	gzip tw_M3_chr${CHR}_2.fastq
#	rm /mnt/quick/julius-temp/mapping/tw_M3_chr${CHR}.bam
#	
#	# re-map them
#	bwa mem -t 9 ~/Documents/mieziai/old/called_M3/consensus_M3_AII-11_chr${CHR}.fa \
#		AII_M3_chr${CHR}_1.fastq.gz AII_M3_chr${CHR}_2.fastq.gz \
#		2> stats/bwa_cons_M3_chr${CHR}_log.txt | samtools fixmate -@2 -m - - | \
#		samtools sort -@4 -m 3500M \
#		-o /mnt/quick/julius-temp/mapping/AII_M3_cons_chr${CHR}.bam
#	bwa mem -t 9 ~/Documents/mieziai/old/called_M3/consensus_M3_AII-11_chr${CHR}.fa \
#		tw_M3_chr${CHR}_1.fastq.gz tw_M3_chr${CHR}_2.fastq.gz \
#		2> stats/bwa_cons_M3_tw_chr${CHR}_log.txt | samtools fixmate -@2 -m - - | \
#		samtools sort -@4 -m 3500M \
#		-o /mnt/quick/julius-temp/mapping/tw_M3_cons_chr${CHR}.bam
#done


# 10. CALL FROM THE RE-MAPS ON CONSENSUS
cd ~/Documents/mieziai/2022/called/
echo "calling on consensus"

#for CHR in {1..7}
#do
#	bcftools mpileup \
#	    -Ou -f ~/Documents/mieziai/old/called_M3/consensus_M3_AII-11_chr${CHR}.fa \
#	    /mnt/quick/julius-temp/mapping/AII_M3_cons_chr${CHR}.bam |\
#	    bcftools call -m -v -Oz -o vcfs/cons/tmp_calledcons_M3_AII_chr${CHR}.vcf.gz
#	bcftools mpileup \
#	    -Ou -f ~/Documents/mieziai/old/called_M3/consensus_M3_AII-11_chr${CHR}.fa \
#	    /mnt/quick/julius-temp/mapping/tw_M3_cons_chr${CHR}.bam |\
#	    bcftools call -m -v -Oz -o vcfs/cons/tmp_calledcons_M3_tw_chr${CHR}.vcf.gz
#	
#	# parse calls, calculate MAF
#	echo 'POS REF ALT DP MQ MQ0F DP4 MAF GT' > callsconsf_AII_M3_chr${CHR}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		vcfs/cons/tmp_calledcons_M3_AII_chr${CHR}.vcf.gz | \
#		awk '{split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' >> callsconsf_AII_M3_chr${CHR}.txt
#	echo 'POS REF ALT DP MQ MQ0F DP4 MAF GT' > callsconsf_tw_M3_chr${CHR}.txt
#	bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#		vcfs/cons/tmp_calledcons_M3_tw_chr${CHR}.vcf.gz | \
#		awk '{split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' >> callsconsf_tw_M3_chr${CHR}.txt
#
#done



# 11. INDIVIDUAL-ONLY COMPARISON
cd ~/Documents/mieziai/old/called_M3/

# extract snps:
bcftools index vcfs/called_M3_tw-12_chr5.vcf.gz
bcftools view --min-af 1 -r 5H:530000000-580000000 -Oz \
	-o annots/called_M3_tw-12_mono.vcf.gz \
	vcfs/called_M3_tw-12_chr5.vcf.gz
# snpEff annotate:
java -jar /mnt/hdd/soft/snpEff/snpEff.jar Hordeum_vulgare \
	annots/called_M3_tw-12_mono.vcf.gz | gzip > annots/annot_M3_tw-12_mono.vcf.gz
# apply some filtering and reformat
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%DP\t%DP4\t%MQ0F\t%MQ\t%ANN\t[%GT]\n" \
	annots/annot_M3_tw-12_mono.vcf.gz | \
       	awk -v FS="\t" 'BEGIN{print "CHR", "POS", "REF", "ALT", "DP", "DP4", "MAF", "MQ0F", "MQ", "ANN"}
	$5>10 && $8>20 && $10=="1\/1"{split($6, d, ","); ref=d[1]+d[2]; alt=d[3]+d[4];
	print substr($1, 1,1), $2, $3, $4, $5, ref+alt, alt/(ref+alt), $7, $8, $9}' > annots/annot_M3_tw-12_mono.txt



exit


## TESTING

# ALTERNATIVE CONSENSUS (only from homozygous gts)
cd ~/Documents/mieziai/2022/mapped/
#bwa mem -t 7 ~/Documents/mieziai/old/called/consensus_AII-11_chr${CHR}_hom.fa \
#	AII_HVgp_nodups_chr${CHR}_1.fastq.gz AII_HVgp_nodups_chr${CHR}_2.fastq.gz \
#	2> bwa_conshom_chr${CHR}_log.txt | samtools fixmate -@2 -m - - | \
#	samtools sort -@4 -m 3500M -o /mnt/quick/julius-temp/mapping/AII_HVgp_conshom_chr${CHR}.bam
#bwa mem -t 7 ~/Documents/mieziai/old/called/consensus_AII-11_chr${CHR}_hom.fa \
#	tw_HVgp_nodups_chr${CHR}_1.fastq.gz tw_HVgp_nodups_chr${CHR}_2.fastq.gz \
#	2> bwa_conshom_tw_chr${CHR}_log.txt | samtools fixmate -@2 -m - - | \
#	samtools sort -@4 -m 3500M -o /mnt/quick/julius-temp/mapping/tw_HVgp_conshom_chr${CHR}.bam
##done

# TESTING

# call from the re-maps on consensus
cd ~/Documents/mieziai/2022/called/
#for CHR in {4..6}
#do
CHR=5
#bcftools mpileup \
#    -Ou -f ~/Documents/mieziai/old/called/consensus_AII-11_chr${CHR}_hom.fa \
#    /mnt/quick/julius-temp/mapping/AII_HVgp_conshom_chr${CHR}.bam |\
#    bcftools call -m -v -Oz -o tmp_calledconshom_AII_chr${CHR}.vcf.gz
#bcftools mpileup \
#    -Ou -f ~/Documents/mieziai/old/called/consensus_AII-11_chr${CHR}_hom.fa \
#    /mnt/quick/julius-temp/mapping/tw_HVgp_conshom_chr${CHR}.bam |\
#    bcftools call -m -v -Oz -o tmp_calledconshom_tw_chr${CHR}.vcf.gz
#
## parse calls
#echo 'POS REF ALT DP MQ MQ0F DP4 GT' > callsconshom_AII_HVgp_chr${CHR}.txt
#bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#	tmp_calledconshom_AII_chr${CHR}.vcf.gz >> callsconshom_AII_HVgp_chr${CHR}.txt
#awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' callsconshom_AII_HVgp_chr${CHR}.txt > callsconshomf_AII_HVgp_chr${CHR}.txt
#
#echo 'POS REF ALT DP MQ MQ0F DP4 GT' > callsconshom_tw_HVgp_chr${CHR}.txt
#bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#	tmp_calledconshom_tw_chr${CHR}.vcf.gz >> callsconshom_tw_HVgp_chr${CHR}.txt
#awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' callsconshom_tw_HVgp_chr${CHR}.txt > callsconshomf_tw_HVgp_chr${CHR}.txt

# pooled consensus
#bcftools mpileup \
#    -Ou -f ~/Documents/mieziai/2022/called/consensus_AII_chr${CHR}.fa \
#    /mnt/quick/julius-temp/mapping/tw_HVgp_conspool_chr${CHR}.bam |\
#    bcftools call -m -v -Oz -o tmp_calledconspool_tw_chr${CHR}.vcf.gz
#
#echo 'POS REF ALT DP MQ MQ0F DP4 GT' > callsconspool_tw_HVgp_chr${CHR}.txt
#bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#	tmp_calledconspool_tw_chr${CHR}.vcf.gz >> callsconspool_tw_HVgp_chr${CHR}.txt
#awk 'NR==1{print $1,$2,$3,$4,$5,$6,"DP4","MAF",$8; next} {split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' callsconspool_tw_HVgp_chr${CHR}.txt > callsconspoolf_tw_HVgp_chr${CHR}.txt

# Duplicate removal check
cd ~/Documents/mieziai/2022/called/
#samtools markdup -@3 -r -s \
#    /mnt/quick/julius-temp/mapping/tw_M3_cons_chr5.bam \
#    /mnt/quick/julius-temp/mapping/tw_M3nodups_cons_chr5.bam &> |
#    /mnt/quick/julius-temp/mapping/dup_stats_tw_M3_cons_chr5.txt
#bcftools mpileup \
#    -Ou -f ~/Documents/mieziai/old/called_M3/consensus_M3_AII-11_chr5.fa \
#    /mnt/quick/julius-temp/mapping/tw_M3nodups_cons_chr5.bam |\
#    bcftools call -m -v -Oz -o vcfs/cons/tmp_calledcons_M3nodups_tw_chr5.vcf.gz
#echo 'POS REF ALT DP MQ MQ0F DP4 MAF GT' > callsconsf_tw_M3nodups_chr5.txt
#bcftools query -f '%POS %REF %ALT %DP %MQ %MQ0F %DP4 [%GT ]\n' \
#	vcfs/cons/tmp_calledcons_M3nodups_tw_chr5.vcf.gz | \
#	awk '{split($7,d,","); rr=d[1]+d[2]; aa=d[3]+d[4]; print $1, $2, $3, $4, $5, $6, rr+aa, aa/(rr+aa), $8}' >> callsconsf_tw_M3nodups_chr5.txt

# MANTA SV calling:
# (separate samples)
cd ~/Documents/mieziai/old/svs/
#MANTA=/mnt/hdd/soft/manta/manta-1.6.0.centos6_x86_64/bin
#$MANTA/configManta.py \
#	--referenceFasta /mnt/quick/julius-temp/Hordeum_vulgare_goldenpromise.GPv1.dna.toplevel.fa \
#	--bam /mnt/quick/julius-temp/mapping/AII-11_HVgp.bam \
#	--runDir ./
