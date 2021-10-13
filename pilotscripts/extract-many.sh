#!/bin/bash
set -e
# Extracts ORFs or other associated regions
# from multiple genes, for REF, WT and tw mutants

# Requires input REF fasta
# curl -O http://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/hordeum_vulgare/dna/Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa.gz
# or equivalent via globus

# Also gff annotations:
# curl -O http://ftp.ensemblgenomes.org/pub/plants/release-51/gff3/hordeum_vulgare/Hordeum_vulgare.IBSC_v2.51.gff3.gz

# And a gene list:
INLIST="auxin_biosynthesis_genes_1.csv"

# and input with all called variants (using those from GATK HaplotypeCaller here).

# extract a gene list from SnpEff database:
# java -jar /mnt/unixarx/julius/soft/snpEff/snpEff.jar dump Hordeum_vulgare -v -txt > snpeff_data_dump.txt
# extract the relevant genes:
# grep -F --color=never -f $INLIST snpeff_data_dump.txt > snpeff_relevant_genes.txt
# grep --color=auto Gene snpeff_relevant_genes.txt > snpeff_onlygenes.txt
# rm snpeff_data_dump.txt

# the output has a list of genes and chromosomes:
# (header is:)
# chr	start	end	strand	type	id	geneName	geneId	numberOfTranscripts	canonicalTranscriptLength	transcriptId	cdsLength	numberOfExons	exonRank	exonSpliceType

# --- MAIN LOOP ---
# for each gene:
i=1
while read chrs starts ends strand type names remainder || [[ -n ${remainder} ]]
do
	echo "--- extracting gene $i : ${names} at chr ${chrs} : ${starts} - ${ends} ---"
	# buffer the upstream direction w/ 1000 bp
	if [ ${strand} == "+1" ]
	then
		starts=$((starts - 1000))
	else
		ends=$((ends + 1000))
	fi

	# extract all variants for that region
	bcftools view -r ${chrs}:${starts}-${ends} ../called/snpeff_11_all.vcf.gz -Oz -o tmp/vars_11_${i}.vcf.gz
	bcftools view -r ${chrs}:${starts}-${ends} ../called/snpeff_12_all.vcf.gz -Oz -o tmp/vars_12_${i}.vcf.gz
	tabix -C tmp/vars_11_${i}.vcf.gz
	tabix -C tmp/vars_12_${i}.vcf.gz
	
	echo "building consensus..."
	# build consensus seqs from ref fasta and variants for each sample
	
	# reference (Morex) - just extract the region
	# and also rename "chr1H" to "1H" to match calls
	samtools faidx ../Hordeum_vulgare.IBSC_v2.dna_sm.toplevel.fa.gz chr${chrs}:${starts}-${ends} | sed 's/chr//'  > extracted/gene${i}_ref.fa
	# WT (ind 11) - compare it w/ the calls
	# tw2 from Auksiniai II (ind 12) - compare it w/ the calls
	cat extracted/gene${i}_ref.fa | bcftools consensus tmp/vars_11_${i}.vcf.gz -s 11 -o extracted/gene${i}_WT.fa
	cat extracted/gene${i}_ref.fa | bcftools consensus tmp/vars_12_${i}.vcf.gz -s 12 -o extracted/gene${i}_tw.fa

	# increment gene counter
	i=$((i + 1))

done < snpeff_onlygenes.txt

# NOTE:
# genes 59-61 are on "Un" chromosome = pieces unassigned to contigs
# so we don't have variant calls for those
# And one of the genes is missing (HORVU6Hr1G0816706) - it is also 1 char longer than all others
