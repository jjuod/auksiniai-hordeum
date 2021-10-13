declare -a genes=("HORVU7Hr1G106280" "HORVU5Hr1G106190" "HORVU2Hr1G076060" "HORVU3Hr1G000700")
declare -a chrs=(7 5 2 3)
declare -a starts=(618942335 623651241 548131773 1737585)
declare -a ends=(618948782 623653262 548134822 1741998)
# get their coords from snpEff

# extract their variants from the called unfiltered and filtered stages
bcftools view -r chr${chrs[0]}H:${starts[0]}-${ends[0]} ../pileup/snpeff_all.vcf.gz -Oz -o gene1_all.vcf.gz
bcftools view -r chr${chrs[0]}H:${starts[0]}-${ends[0]} ../pileup/snpeff_filtered.vcf.gz -Oz -o gene1_filt.vcf.gz
bcftools view -r chr${chrs[0]}H:${starts[0]}-${ends[0]} ../pileup/snpeff_gt_filtered.vcf.gz -Oz -o gene1_gt_filt.vcf.gz

bcftools view -r chr${chrs[1]}H:${starts[1]}-${ends[1]} ../pileup/snpeff_all.vcf.gz -Oz -o gene2_all.vcf.gz
bcftools view -r chr${chrs[1]}H:${starts[1]}-${ends[1]} ../pileup/snpeff_filtered.vcf.gz -Oz -o gene2_filt.vcf.gz
bcftools view -r chr${chrs[1]}H:${starts[1]}-${ends[1]} ../pileup/snpeff_gt_filtered.vcf.gz -Oz -o gene2_gt_filt.vcf.gz

bcftools view -r chr${chrs[2]}H:${starts[2]}-${ends[2]} ../pileup/snpeff_all.vcf.gz -Oz -o gene3_all.vcf.gz
bcftools view -r chr${chrs[2]}H:${starts[2]}-${ends[2]} ../pileup/snpeff_filtered.vcf.gz -Oz -o gene3_filt.vcf.gz
bcftools view -r chr${chrs[2]}H:${starts[2]}-${ends[2]} ../pileup/snpeff_gt_filtered.vcf.gz -Oz -o gene3_gt_filt.vcf.gz

bcftools view -r chr${chrs[3]}H:${starts[3]}-${ends[3]} ../pileup/snpeff_all.vcf.gz -Oz -o gene4_all.vcf.gz
bcftools view -r chr${chrs[3]}H:${starts[3]}-${ends[3]} ../pileup/snpeff_filtered.vcf.gz -Oz -o gene4_filt.vcf.gz
bcftools view -r chr${chrs[3]}H:${starts[3]}-${ends[3]} ../pileup/snpeff_gt_filtered.vcf.gz -Oz -o gene4_gt_filt.vcf.gz

# build consensus seqs from ref fasta and variants for each indiv, using only the FULL files
for g in 0 1 2 3
do
    echo "region ${chrs[$g]}: ${starts[$g]} - ${ends[$g]}"
    gi=$(($g + 1))
    echo "gi = $gi"
    tabix -C gene${gi}_all.vcf.gz
    # reference (Morex)
    samtools faidx ../Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa chr${chrs[$g]}H:${starts[$g]}-${ends[$g]} > gene${gi}_ref.fa
    # WT and tw from Auksiniai II
    samtools faidx ../Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa chr${chrs[$g]}H:${starts[$g]}-${ends[$g]} | bcftools consensus gene${gi}_all.vcf.gz -s 11 -o gene${gi}_11.fa
    samtools faidx ../Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa chr${chrs[$g]}H:${starts[$g]}-${ends[$g]} | bcftools consensus gene${gi}_all.vcf.gz -s 12 -o gene${gi}_12.fa
done

# extract info from snpEff tables
awk '$0~/Exons/{p=1;next} $0~/CDS/{p=2} $0~/UTR/{p=2} p==1{split($1, a, ":"); split(a[2], b, "-"); print("chr" a[1], b[1], b[2])}' extracted_g*.vcf > exons.txt
