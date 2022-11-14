library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/Documents/mieziai/calls/genes/")

# Extracting the target chr 5 region SNPs
# from morex v3 annotation

# snpEff annotations for the tw SNPs
ann = fread("annot_target_M3_tw_mono.txt")
# already filtered by GT, MQ, DP
# filter strictly by MAF:
ann = filter(ann, MAF>0.9, MQ>30, DP>=15, DP<90)

ann = mutate(ann, conseq = ifelse(grepl("HIGH", ANN), "high",
                                  ifelse(grepl("MODERATE", ANN), "mod", "low"))) %>%
  mutate(conseq = factor(conseq, levels=c("high", "mod", "low")))
table(ann$conseq)

# AII pool SNPs
callsA = fread("../callsf_AII_M3_chr5.txt")
callsA = filter(callsA, POS>529e6, POS<575e6)
nrow(callsA)  # 565 k

# Filtering for mapping quality:
callsA = filter(callsA, MQ>30)
nrow(callsA)  # 496 k

# Filtering by depth:
callsA = filter(callsA, DP>=15, DP<90)
nrow(callsA)  # 451 k

callsA$id = paste(callsA$POS, callsA$REF, callsA$ALT, sep="_")

ann = mutate(ann, in_contr = paste(POS, REF, ALT, sep="_") %in% callsA$id)
table(ann$in_contr, ann$conseq)


# for fun: load GP hit snps for comparison
annot_tw = fread("../../annots/annot_tw_mono.txt")
annot_tw = filter(annot_tw, POS>490e6, POS<530e6, CHR==5)
nrow(annot_tw)

annot_tw = filter(annot_tw, MAF>0.9, MQ>30, DP>=15, DP<90)
annot_tw = mutate(annot_tw, conseq = ifelse(grepl("HIGH", ANN), "high",
                                  ifelse(grepl("MODERATE", ANN), "mod", "low"))) %>%
  mutate(conseq = factor(conseq, levels=c("high", "mod", "low")))
nrow(annot_tw)
table(annot_tw$conseq)

# interesting:
# in GP:
# high  mod      low 
# 76   1365   118509 
# in M3:
# high  mod      low 
# 79   1843   204214 

callsA_GP = fread("../callsf_AII_HVgp_nodups_chr5.txt")
callsA_GP = filter(callsA_GP, POS>490e6, POS<530e6)
callsA_GP = filter(callsA_GP, MQ>30, DP>=15, DP<90)
nrow(callsA_GP)  # 404 k

callsA_GP$id = paste(callsA_GP$POS, callsA_GP$REF, callsA_GP$ALT, sep="_")

annot_tw = mutate(annot_tw, in_contr = paste(POS, REF, ALT, sep="_") %in% callsA_GP$id)
table(annot_tw$in_contr, annot_tw$conseq)


# Formatting nicely:
ann2 = filter(ann, !in_contr)
ann2 = separate_rows(ann2, ANN, sep=",") %>%
  separate(ANN,
           c("Allele", "Effect", "Impact", "GeneName", "GeneID", "FeatureType", "FeatureID", 
             "TranscrBiotype", "ExonRank", "HVGS.c", "HGVS.p", "pos.cDNA", "pos.CDS", "pos.PROT",
             "Distance", "Errors"), sep="\\|")
ann2 = select(ann2, -one_of(c("GeneID", "FeatureType", "FeatureID", "TranscrBiotype",
                                      "Distance", "conseq", "Allele", "in_contr")))

ann_hi = ann2 %>% filter(Impact=="HIGH")
ann_mod = ann2 %>% filter(Impact=="MODERATE")

nrow(ann_hi)  # 29
nrow(ann_mod)  # 391

# attach GO terms and save
biomart = fread("../../refs/biomart_go_annots_M3.txt")
biomart = unite(biomart, "Pfam.pos", "Pfam start":"Pfam end", sep="-") %>%
  unite("Pfam.domains", "Pfam ID":"Pfam.pos", sep=":")
biomart = group_by(biomart, `Gene stable ID`) %>%
  summarize(GO.term.name=paste(unique(`GO term name`), collapse="; "),
            Pfam.domains=paste(unique(`Pfam.domains`), collapse="; "))

ann_hi_go = biomart %>%
  inner_join(ann_hi, ., by=c("GeneName"="Gene stable ID"))

ann_mod_go = biomart %>%
  inner_join(ann_mod, ., by=c("GeneName"="Gene stable ID"))

print.data.frame(ann_hi_go)
write.table(ann_hi_go, "~/Documents/mieziai/calls/final_tw_M3_chr5_hi.csv", sep="\t", quote=F, row.names=F)
write.table(ann_mod_go, "~/Documents/mieziai/calls/final_tw_M3_chr5_mod.csv", sep="\t", quote=F, row.names=F)

