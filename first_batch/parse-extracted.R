# HORVU7Hr1G106280	4341841	chr7H	618948243

# HORVU5Hr1G106190	chr5H	623651732	TCC	p.Pro129fs,p.Pro129fs,p.Pro79fs,p.Pro80fs,p.Pro93fs,p.Pro14fs		tryptophan synthase activity

# HORVU2Hr1G076060	chr2H	548134798	C	p.Gln9Glu			amino acid transport

# HORVU3Hr1G000700	chr3H	1737837	G	p.Ala49fs,p.Ala84fs	WAT1-related protein

options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
setwd("~/Documents/mieziai/")

seqs = data.frame()

# ref fastas
for(g in 1:4){
  ref = readLines(paste0("extractedgenes/gene", g, "_ref.fa"))
  newseq = data.frame(gene=g, header=ref[1], ind="ref", s=paste(ref[-1], collapse=""))
  seqs = bind_rows(seqs, newseq)
}

# variant vcfs
vars = data.frame()
effs = data.frame()
for(g in 1:4){
  allvars = read.table(gzfile(paste0("extractedgenes/gene", g, "_all.vcf.gz")))
  allvars = separate(allvars, V10, c("GT11", "PL11"), sep=":")
  allvars = separate(allvars, V11, c("GT12", "PL12"), sep=":")
  # so just separate out depth, depth per allele, and annotation
  allvars = extract(allvars, V8, c("DP", "DP4", "ANN"),
                  regex=c("DP=([0-9]*).*DP4=([0-9,]*).*ANN=(.*)"),
                  convert=T)
  allvars$gene = g
  effs = bind_rows(effs, allvars[,c(1:5, 10)])
  vars = select(allvars, -one_of("V9", "ANN", "PL11", "PL12", "V3", "V7")) %>%
                  bind_rows(vars, .)
}
effs = effs[,-3]
effs$worst = ifelse(grepl("HIGH", effs$ANN), "HIGH", ifelse(grepl("MODERATE", effs$ANN), "MODERATE", "LOW"))

# filtered and genotype-filtered (00 WT, 11 tw) vcfs
filt = data.frame()
gtfilt = data.frame()
for(g in 1:4){
  filtvars = read.table(gzfile(paste0("extractedgenes/gene", g, "_filt.vcf.gz")))
  filt = bind_rows(filt, filtvars[,c(1, 2, 4, 5)])
  filtvars = read.table(gzfile(paste0("extractedgenes/gene", g, "_gt_filt.vcf.gz")))
  gtfilt = bind_rows(gtfilt, filtvars[,c(1, 2, 4, 5)])
}

filt$passfilt = T
gtfilt$passgt = T
vars = left_join(vars, filt, by=c("V1", "V2", "V4", "V5"))
vars = left_join(vars, gtfilt, by=c("V1", "V2", "V4", "V5"))
vars$passfilt = !is.na(vars$passfilt)
vars$passgt = !is.na(vars$passgt)
vars$pass = ifelse(vars$passfilt, ifelse(vars$passgt, "PASS", "GTFAIL"), "FILTFAIL")
table(vars$passfilt)
table(vars$passgt)

exons = read.table("extractedgenes/exons.txt")
exons = unique(exons)
nrow(exons)
exons$id = seq_along(exons$V1)/nrow(exons)

# PLOT
filter(vars, gene==1) %>%
  gather(key="ind", value="GT", GT11:GT12) %>%
  ggplot() + geom_point(aes(x=V2, col=pass, y=ind, size=GT), pch="|")+
  geom_segment(data=filter(exons, V1==vars$V1[1]), aes(x=V2, xend=V3, y="ex", yend="ex")) +
  scale_color_discrete(limits=c("FILTFAIL", "PASS", "GTFAIL")) +
  scale_size_manual(limits=c("./.", "0/0", "0/1", "1/1"), values = c(2, 0.5, 3, 6)) +
  theme_bw()

# PLOT only pass variants
filter(vars, gene==1, passgt, passfilt) %>%
  left_join(effs, by=c("V1", "V2", "V4", "V5")) %>% View
filter(vars, gene==1, passgt, passfilt) %>%
  left_join(effs, by=c("V1", "V2", "V4", "V5")) %>%
  ggplot() + geom_point(aes(x=V2, col=worst), y=1, size=7, pch="|") +
  geom_segment(data=filter(exons, V1==vars$V1[1]), aes(x=V2, xend=V3, y=id, yend=id)) +
  theme_bw() + coord_cartesian(ylim=c(0, 2))
