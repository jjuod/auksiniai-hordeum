set -e

M3REF=/mnt/quick/julius-temp/Hordeum_vulgare.MorexV3.dna_sm.fa 

# 5. map onto ref genome

# . remove duplicates, while doing appropriate sorts
# NOTE: not done, the duplicate number was small on these.

# DYSGU SV calling:
cd ~/Documents/mieziai/2024/svs/
dysgu run ${M3REF} -p14 ../tmp/ \
	~/Documents/mieziai/old/mapped_M3/AII-11_M3.bam \
	-x > M3/dysgu_AII-11.vcf
dysgu run ${M3REF} -p14 ../tmp/ \
	~/Documents/mieziai/old/mapped_M3/tw-12_M3.bam \
	-x > M3/dysgu_tw-12.vcf
# compare w/ wildtype and then filter on quality
# (also tried the other way around, got similar amounts of svs)
dysgu filter --normal-vcf M3/dysgu_AII-11.vcf \
	M3/dysgu_tw-12.vcf \
	~/Documents/mieziai/old/mapped_M3/AII-11_M3.bam > M3/dysgu_tw-12_unique.vcf
dysgu filter --min-prob 0.2 --support-fraction 0.15 \
	--pass-prob 0.3 \
	M3/dysgu_tw-12_unique.vcf > M3/dysgu_tw-12_unique_filteredafter.vcf


# DYSGU on pools:
cd ~/Documents/mieziai/2024/svs/
#dysgu run ${M3REF} -p14 ../tmp/ \
#	~/Documents/mieziai/2022/mapped/AII_M3.bam \
#	-x > M3/dysgu_AII.vcf
dysgu run ${M3REF} -p14 ../tmp/ \
	~/Documents/mieziai/2022/mapped/tw_M3.bam \
	-x > M3/dysgu_tw.vcf
# compare w/ wildtype and then filter on quality
dysgu filter --normal-vcf M3/dysgu_AII.vcf \
	M3/dysgu_tw.vcf \
	~/Documents/mieziai/2022/mapped/AII_M3.bam > M3/dysgu_tw_unique.vcf
dysgu filter --min-prob 0.2 --support-fraction 0.15 \
	--pass-prob 0.3 \
	M3/dysgu_tw_unique.vcf > M3/dysgu_tw_unique_filteredafter.vcf

echo "SV calling via DYSGU complete."

# MANTA SV calling:
# (separate samples)
cd ~/Documents/mieziai/old/svs/

# generate configs
MANTA=/mnt/hdd/soft/manta/manta-1.6.0.centos6_x86_64/bin
#$MANTA/configManta.py \
#	--referenceFasta ${M3REF} \
#	--bam ~/Documents/mieziai/old/mapped_M3/AII-11_M3.bam \
#	--runDir ./AII-11_M3/
#
#$MANTA/configManta.py \
#	--referenceFasta ${M3REF} \
#	--bam ~/Documents/mieziai/old/mapped_M3/tw-12_M3.bam \
#	--runDir ./tw-12_M3/
$MANTA/configManta.py \
	--normalBam ~/Documents/mieziai/old/mapped_M3/AII-11_M3.bam \
	--tumorBam ~/Documents/mieziai/old/mapped_M3/tw-12_M3.bam \
	--referenceFasta ${M3REF} \
	--runDir ./tw-12_specific_M3/

#./AII-11_M3/runWorkflow.py
#./tw-12_M3/runWorkflow.py
./tw-12_specific_M3/runWorkflow.py
#
# # manta breaks at the last step (making indices), 
# # run tabix -C on the files if need indices
#echo "SV calling via MANTA complete."

exit


# ------------ NOT RUN ---------------
# Duplicate removal check
#cd ~/Documents/mieziai/2022/called/
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

