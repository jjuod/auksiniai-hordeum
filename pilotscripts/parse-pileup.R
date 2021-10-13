options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)

# Script for parsing samtools mpileup + bcftools call style results.

# read in
calls = read.table(gzfile("~/Documents/mieziai/pileup/snpeff_gt_filtered.vcf.gz"))

# 16 M known variants, filtered for the same positions
knownvars = read.table(gzfile("~/Documents/mieziai/pileup/ensembl-vars_gt_filtered.vcf.gz"), sep="\t")
#knownvars = filter(knownvars, V2<max(calls$V2))

# test file w/o PCR duplicates
nodup = read.table(gzfile("~/Documents/mieziai/pileup/nodup_called.vcf.gz"))
potdup = filter(calls, V2<10000000)
mean(potdup$known | potdup$knownRC)
nodup = semi_join(potdup, nodup, by=c("V1", "V2", "V4", "V5"))
mean(nodup$known | nodup$knownRC)

# split columns - does not exactly work because some columns are optional
#calls = extract(calls, V8, c("DP", "VDB", "SGB", "RPB", "MQB",
#                             "MQSB", "BQB", "MQ0F", "ICB", "HOB",
#                             "AC", "AN", "DP4", "MQ", "ANN"),
#                regex="DP=(.*);VDB=(.*);SGB=(.*);RPB=(.*);MQB=(.*);MQSB=(.*);BQB=(.*);MQ0F=(.*);ICB=(.*);HOB=(.*);AC=(.*);AN=(.*);DP4=(.*);MQ=(.*);ANN=(.*)")

# so just separate out depth, depth per allele, and annotation
calls = extract(calls, V8, c("DP", "DP4", "ANN"),
                regex=c("DP=([0-9]*).*DP4=([0-9,]*).*ANN=(.*)"),
                convert=T)

# split geno columns into GT (hardcalls) and PL (likelihoods)
calls = separate(calls, V10, c("GT11", "PL11"), sep=":")
calls = separate(calls, V11, c("GT12", "PL12"), sep=":")


## MERGE WITH KNOWN VARIATION

# produce reverse complement cause many strands are flipped:
knownvarsRC = knownvars
knownvarsRC$V4[knownvars$V4=="A"] = "T"
knownvarsRC$V4[knownvars$V4=="T"] = "A"
knownvarsRC$V4[knownvars$V4=="C"] = "G"
knownvarsRC$V4[knownvars$V4=="G"] = "C"

knownvarsRC$V5[knownvars$V5=="A"] = "T"
knownvarsRC$V5[knownvars$V5=="T"] = "A"
knownvarsRC$V5[knownvars$V5=="C"] = "G"
knownvarsRC$V5[knownvars$V5=="G"] = "C"

# collapse imported df
knownvars = paste(knownvars$V1, knownvars$V2, knownvars$V4, knownvars$V5, sep="_")
knownvarsRC = paste(knownvarsRC$V1, knownvarsRC$V2, knownvarsRC$V4, knownvarsRC$V5, sep="_")

# attach flags to main file
calls = mutate(calls, known = paste(V1, V2, V4, V5, sep="_") %in% knownvars,
       knownRC = paste(V1, V2, V4, V5, sep="_") %in% knownvarsRC)
table(calls$known, calls$knownRC)

length(knownvars) # 265 k known
sum(calls$known) + sum(calls$knownRC) # 94 k found in our thing
nrow(calls) # out of 500 k calls
mean(calls$known) + mean(calls$knownRC) # ~20 % of variants already known, even if we merge conservatively (no indels...)

### CHECK GENOTYPES ACROSS STRAINS

# get table of genotypes by line
byline = group_by(calls, GT11, GT12) %>%
  summarize(n=n(), fr=n/nrow(calls), mdp=mean(DP), known=mean(known)+mean(knownRC)) %>%
  arrange(desc(n))
byline %>%
  View()

group_by(nodup, GT11, GT12) %>%
  summarize(n=n(), fr=n/nrow(calls), mdp=mean(DP), known=mean(known)+mean(knownRC)) %>%
  arrange(desc(n)) %>% View()

# most of the known variants detected here are 1/1-1/1 calls
# or: 28 % of our 1/1-1/1 calls are already known
# compare w/ 4 % of our 0/1-0/1 calls, or 5 % of 1/1-./. calls
# BUT: 33 % of our 0/0-1/1 calls known??

alleles = group_by(calls, V4) %>%
  summarize(n = n(), size = nchar(as.character(V4[1])))
alleles %>%
  arrange(desc(n)) %>%
  View()

### FILTER BY DEPTH
group_by(calls, DP=round(DP)) %>%
  summarize(known = mean(known | knownRC), n=n()) %>%
  filter(n>10) %>%
  ggplot(aes(x=DP, y=known)) + geom_line()

# peak of known variation seen at DP 50
group_by(calls, DP>50) %>%
  summarize(known = mean(known | knownRC), n=n())
# everything with DP<=10 can be safely filtered out I think
group_by(calls, DP>10) %>%
  summarize(known = mean(known | knownRC), n=n())

callsF = filter(calls, DP>10)
group_by(callsF, GT11, GT12) %>%
  summarize(n=n(), fr=n/nrow(callsF), mdp=mean(DP), known=mean(known)+mean(knownRC)) %>%
  arrange(desc(n)) %>% View()

callsF = filter(calls, GT11=="0/0", GT12=="1/1", !known, !knownRC)
table(grepl("HIGH", callsF$ANN))
table(grepl("MODERATE", callsF$ANN))

filter(callsF, grepl("HIGH", ANN))

# split by gene and drop empty rows at line ending
#calls = separate_rows(calls, ANN, sep=",") %>%
  #filter(ANN!="")

nrow(calls)


### FINAL FILTERING
# read in
calls = read.table(gzfile("~/Documents/mieziai/pileup/snpeff_gt_filtered.vcf.gz"))

# known variants, filtered for the same positions
knownvars = read.table(gzfile("~/Documents/mieziai/pileup/ensembl-vars_gt_filtered.vcf.gz"), sep="\t")

# separate out depth, fraction of mapQ=0 reads, depth per allele, rms mapQ, and annotation
calls = extract(calls, V8, c("DP", "MQ0F", "DP4", "MQ", "ANN"),
                regex=c("DP=([0-9]*).*MQ0F=([^;]*).*DP4=([0-9,]*).*MQ=([^;]*).*ANN=(.*)"),
                convert=T)

# split geno columns into GT (hardcalls) and PL (likelihoods)
calls = separate(calls, V10, c("GT11", "PL11"), sep=":")
calls = separate(calls, V11, c("GT12", "PL12"), sep=":")

# drop ID and INFO (V3 and V7)
calls = select(calls, -one_of(c("V3", "V7", "V9")))

# filter mutations with non-"Modifier" effects, i.e. at least somewhat deleterious
callsDel = filter(calls, grepl("LOW|MODERATE|HIGH", ANN))
sum(grepl("LOW", callsDel$ANN))
sum(grepl("MODERATE", callsDel$ANN))
sum(grepl("HIGH", callsDel$ANN)) 
nrow(calls)

# changes due to filtering - less deleterious mutations
sum(grepl("HIGH", callsDel$ANN)) / nrow(calls) * 100 # 0.057, decrease from 0.155%
sum(!grepl("HIGH|LOW", callsDel$ANN)) / nrow(calls) * 100 # 0.413, decrease from 0.839%
sum(!grepl("HIGH|MODERATE", callsDel$ANN)) / nrow(calls) * 100 # 0.486, decrease from 0.937%
nrow(callsDel) / nrow(calls) * 100 # 1.004

sum(grepl("missense_variant", callsDel$ANN)) / nrow(calls) * 100 # 0.445, decrease from 0.799%
group_by(calls, V4, V5) %>%
  summarize(n=n()) %>%
  filter(n>10000) %>%
  mutate(n = n/nrow(calls)*100)
# transitions:
12.2+16.4+16.4+12.2  
# transversions:
4.31+3.27+5.58+5.30+5.30+5.59+3.24+4.28

# known variation fraction:
# produce reverse complement cause many strands are flipped:
knownvarsRC = knownvars
knownvarsRC$V4[knownvars$V4=="A"] = "T"
knownvarsRC$V4[knownvars$V4=="T"] = "A"
knownvarsRC$V4[knownvars$V4=="C"] = "G"
knownvarsRC$V4[knownvars$V4=="G"] = "C"

knownvarsRC$V5[knownvars$V5=="A"] = "T"
knownvarsRC$V5[knownvars$V5=="T"] = "A"
knownvarsRC$V5[knownvars$V5=="C"] = "G"
knownvarsRC$V5[knownvars$V5=="G"] = "C"

# collapse imported df
knownvars = paste(knownvars$V1, knownvars$V2, knownvars$V4, knownvars$V5, sep="_")
knownvarsRC = paste(knownvarsRC$V1, knownvarsRC$V2, knownvarsRC$V4, knownvarsRC$V5, sep="_")

# attach flags to main file
callsDel = mutate(callsDel, known = paste(V1, V2, V4, V5, sep="_") %in% knownvars,
               knownRC = paste(V1, V2, V4, V5, sep="_") %in% knownvarsRC)
sum(callsDel$knownRC)
sum(callsDel$known)

# Remove known variation: (~33%)
callsDel = filter(calls, grepl("LOW|MODERATE|HIGH", ANN))
callsDel = filter(callsDel, !paste(V1, V2, V4, V5, sep="_") %in% knownvars,
                  !paste(V1, V2, V4, V5, sep="_") %in% knownvarsRC)

callsHi = filter(callsDel, grepl("HIGH", ANN))
callsHiMod = filter(callsDel, grepl("HIGH|MODERATE", ANN))
nrow(callsHi) # 653
nrow(callsHiMod) # 3477

qplot(callsHi$DP)
qplot(callsHi$MQ0F)

# export positions for confirming with HaplotypeCaller
mutate(callsHiMod, V1 = substr(V1, 4, 6))[,c("V1", "V2")] %>%
  write.table(., "~/Documents/mieziai/pileup/positions_mod-hi.tsv",
            sep="\t", col.names=F, row.names=F, quote=F)

rm(calls)

# HaplotypeCaller annots
hcHiMod1 = read.table(gzfile("~/Documents/mieziai/called/overlHC_11_mod-hi.vcf.gz"))
hcHiMod2 = read.table(gzfile("~/Documents/mieziai/called/overlHC_12_mod-hi.vcf.gz"))

# MOST of the annotations match:
table(grepl("HIGH|MODERATE", hcHiMod2$V8))
hcHiMod1 = separate(hcHiMod1, V10, c("GT11hc", "AD11hc", "DP11hc", "GQ11hc", "PL11hc"), sep=":")
hcHiMod2 = separate(hcHiMod2, V10, c("GT12hc", "AD12hc", "DP12hc", "GQ12hc", "PL12hc"), sep=":")
table(hcHiMod1$GT11hc)
table(hcHiMod2$GT12hc)
nrow(hcHiMod2)

## GOOD SNPS SHOULD:
# have confirmed 12-1/1 in both:
goodHcSnps = filter(hcHiMod2,GT12hc=="1/1")
nrow(goodHcSnps) # 3080

# have no noticeable detection in 11:
goodHcSnps = anti_join(goodHcSnps, hcHiMod1, by=c("V1", "V2", "V4", "V5"))
nrow(goodHcSnps) # 3062

# bad annotations (these should differ only if the called alleles differed slightly)
goodHcSnps = filter(goodHcSnps, grepl("HIGH|MODERATE", V8))
nrow(goodHcSnps) # 3036

goodHcSnps$V1 = paste0("chr", goodHcSnps$V1)
goodHcSnps = goodHcSnps[,c("V1", "V2", "V4", "V5", "DP12hc")]

goodHcSnps = unique(goodHcSnps)
nrow(goodHcSnps) # 3021. shrug

# filter:
callsHiMod$matched = FALSE

for(r in 1:nrow(callsHiMod)){
  # merge on exact chr and +- fuzzy boundary
  matches = filter(goodHcSnps, V1==callsHiMod$V1[r],
         abs(V2-callsHiMod$V2[r]) < max(nchar(V4), nchar(V5)))

  if (nrow(matches)>0){
    callsHiMod$matched[r] = TRUE
  }
}
table(callsHiMod$matched)


### FINAL DATASETS:
callsHiModHC = filter(callsHiMod, matched) %>%
  select(-one_of("matched"))
nrow(callsHiModHC) # 3111

callsHiHC = filter(callsHiModHC, grepl("HIGH", ANN))
nrow(callsHiHC) # 572


## another version with cleaner ANN field, separated per gene
genes = separate_rows(callsHiHC, ANN, sep=",")
genes = genes[,c("V1", "V2", "ANN")]
genes = separate(genes, ANN, c("Allele", "Ann", "Impact", "GeneName", "GeneID", "FeatureType", "FeatureID", "Coding",
                               "ExonIntronRank", "ChangeDNA", "ChangeProt", "PosInTranscript",
                               "PosInCDS", "AAPos", "Distance", "Warnings"), sep="\\|", extra="merge")

genesMod = separate_rows(callsHiModHC, ANN, sep=",")
genesMod = genesMod[,c("V1", "V2", "ANN")]
genesMod = separate(genesMod, ANN, c("Allele", "Ann", "Impact", "GeneName", "GeneID", "FeatureType", "FeatureID", "Coding",
                               "ExonIntronRank", "ChangeDNA", "ChangeProt", "PosInTranscript",
                               "PosInCDS", "AAPos", "Distance", "Warnings"), sep="\\|", extra="merge")

## GO
go = read.table("~/Documents/mieziai/biomart-hv-go.tsv", sep="\t", quote = "", h=T)
go = unique(go)
go = filter(go, GO.term.name!="")

genes_bygo = filter(genes, Impact=="HIGH") %>% 
  group_by(GeneID, GeneName, V1, V2, Allele) %>%
  summarize(Mutations = paste(ChangeProt, collapse=",")) %>%
  left_join(go, by=c("GeneID"="Gene.stable.ID"))

genesMod_bygo = filter(genesMod, Impact=="HIGH" | Impact=="MODERATE") %>% 
  group_by(GeneID, GeneName, V1, V2, Allele) %>%
  summarize(Mutations = paste(ChangeProt, collapse=",")) %>%
  left_join(go, by=c("GeneID"="Gene.stable.ID"))

# SAVE OUTPUTS
colnames(callsHiHC)[1:5] = c("CHR", "POS", "REF", "ALT", "QUAL")
colnames(callsHiModHC)[1:5] = c("CHR", "POS", "REF", "ALT", "QUAL")
write.table(callsHiHC, "~/Documents/mieziai/finalmuts/byvariant-hi.csv", sep="\t", quote=F, row.names=F)
write.table(callsHiModHC, "~/Documents/mieziai/finalmuts/byvariant-himod.csv", sep="\t", quote=F, row.names=F)

colnames(genes)[1:2] = c("CHR", "POS")
colnames(genesMod)[1:2] = c("CHR", "POS")
write.table(genes, "~/Documents/mieziai/finalmuts/bygene-hi.csv", sep="\t", quote=F, row.names=F)
write.table(genesMod, "~/Documents/mieziai/finalmuts/bygene-himod.csv", sep="\t", quote=F, row.names=F)

colnames(genes_bygo)[3:4] = c("CHR", "POS")
colnames(genesMod_bygo)[3:4] = c("CHR", "POS")
write.table(genes_bygo, "~/Documents/mieziai/finalmuts/bygo-hi.csv", sep="\t", quote=F, row.names=F)
write.table(genesMod_bygo, "~/Documents/mieziai/finalmuts/bygo-himod.csv", sep="\t", quote=F, row.names=F)

# barley chr length:
roh = read.table("~/Documents/mieziai/mapped/pileup/roh.txt.gz")
rohwt = read.table("~/Documents/mieziai/mapped/pileup/roh_wt.txt.gz")

ggplot(roh) + geom_segment(aes(x=V4, xend=V5, y=V3, yend=V3),col="red")
roh2 = group_by(roh, V3) %>% summarize(l=sum(V6))
roh2wt = group_by(rohwt, V3) %>% summarize(l=sum(V6))
roh2$lwt = roh2wt$l
roh2$normal = c(558535000, 768075000, 699711000, 647060000, 670030000, 583380000, 657224000)
# c(558535000, 768075000, 699711000, 647060000, 670030000, 583380000, 657224000)
ggplot(roh) +
  facet_wrap(~V3) +
  geom_point(aes(x=V4, y=V6),col="red") +
  geom_point(aes(x=V4, y=V6),col="blue", data=rohwt)

ggplot(filter(roh, V3=="chr6H", V4>1.5e8, V4<2e8)) +
  geom_point(aes(x=V4, y=V6),col="red")


# gromm for variant calling: case/control, valgo .bam
# reiktu palyginti, kokius snpus randa gatk
# gorilla- go term enrichment