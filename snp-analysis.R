library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
setwd("~/Documents/mieziai/calls/")

# VCF fields:
# ID=MQ0F,Type=Float,Description="Fraction of MQ0 reads (smaller is better)"
# ID=MQ,Type=Integer,Description="Average mapping quality"
# ID=DP,Number=1,Type=Integer,Description="Raw read depth"
# ID=DP4,Number=4,Description="Number of high-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases

# Based on diagnostics.Rmd:
# - use GoldenPromise ref
# - duplicates can be removed
# - can use joint calling, but need separate depths anyway

# load genetic map
# from  IBSC (2016-07-12): High-resolution GBS map of the Morex x Barke RIL population. DOI:10.5447/ipk/2016/29
map = read.table("../refs/barley_morex_x_barke_high_resolution_gbs_map.tsv", h=T)
map = map[,2:4]
map$chr = as.numeric(substr(map$chr, 4, 4))

calls_out = tibble()
calls_out_cm = tibble()
for(chr in 1:7){
  print(chr)
  callsA = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                     "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                    .id="pool")
  callsE = bind_rows("AII"=fread(paste0("exons_AII_HVgp_nodups_chr",chr,".txt")),
                    "tw"=fread(paste0("exons_tw_HVgp_nodups_chr",chr,".txt")),
                    .id="pool")
  calls = bind_rows("all"=callsA, "exon"=callsE, .id="type")
  
  # Filtering for mapping quality:
  calls = filter(calls, MQ>30)
  
  # Filtering by depth:
  calls = filter(calls, DP>=15, DP<90)  # about 10 % loss
  
  calls = calls[,c("pool", "type", "POS", "REF", "ALT", "DP", "MAF")]
  
  # Filtering only those present in both samples:
  callsBOTH = inner_join(filter(calls, pool=="AII", type=="all"),
                     filter(calls, pool=="tw", type=="all"),
                     by=c("POS", "REF", "ALT"),
                     suffix=c(".AII", ".tw")) %>%
    select(-one_of(c("pool.AII", "pool.tw")))
  # This is important! otherwise depth etc. can really differ
  # depending on which markers are included

  # Keep only unique ones for either pool:
  callsUNIQ = anti_join(filter(calls, type=="all"),
                        callsBOTH, by=c("POS", "REF", "ALT"))

  # Filtering to remove roughly monomorphic markers:
  callsBOTH = filter(callsBOTH, MAF.AII<0.8 | MAF.tw<0.8)

  callsBOTH = pivot_longer(callsBOTH, c(DP.AII, DP.tw, MAF.AII, MAF.tw),
               names_to=c("stat", "pool"), names_sep="\\.", values_to="val") %>%
    spread("stat", "val")

  callsBOTH$type="shared"
  callsUNIQ$type="unique"
  calls = bind_rows(calls, callsBOTH, callsUNIQ)
  
  calls = mutate(calls, w = POS %/% 2e6 * 2)  # in 2 Mbp windows
  
  # merge in the genetic coords
  map_chr = map[map$chr==chr,]
  # basic rescaling, to at least match the length of GP chroms:
  range_old = range(map_chr$pos)
  range_new = range(calls$POS)
  map_chr$pos = (map_chr$pos-range_old[1]) / diff(range_old) * diff(range_new) + range_new[1]
  calls$cM = approx(map_chr$pos, map_chr$gbs_cM, calls$POS, rule=2)$y
  calls = mutate(calls, w_cM = cM %/% 1 * 1)  # in 1 cM windows
  
  calls_sum = calls %>%
    group_by(w, pool, type) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))
  calls_sum_cm = calls %>%
    group_by(w_cM, pool, type) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))

  calls_out = bind_rows(calls_out, calls_sum)
  calls_out_cm = bind_rows(calls_out_cm, calls_sum_cm)
}

# Physical map:
# remember that depth changes basically only b/c of selection
# on called positions, and is quite flat overall
calls_out %>%
  ggplot(aes(x=w, y=medDP)) +
  geom_line(aes(col=pool), alpha=0.7) +
  facet_grid(type~chr, scales="free_x") +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out, aes(x=w, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_grid(type~chr, scales="free_x") +
  theme_bw() + theme(legend.position="bottom")

# Genetic map:
ggplot(calls_out_cm, aes(x=w_cM, y=meanMAF, col=pool)) +
  geom_line(alpha=0.4) +
  facet_grid(type~chr, scales="free_x") +
  geom_smooth(se=F) + 
  theme_bw() + theme(legend.position="bottom")

# Frac AF>0.9, instead of mean AAF, cleans up the noise quite well
# as many positions have moderate AFs:
ggplot(calls_out_cm, aes(x=w_cM, y=mean11, col=pool)) +
  geom_line(alpha=0.4) +
  facet_grid(type~chr, scales="free_x") +
  geom_smooth(se=F) + 
  theme_bw() + theme(legend.position="bottom")

# --------


# MAFs and median depth of called variants
# over 2 Mbp physical windows
ggplot(calls_out, aes(x=w, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggsave("~/Documents/mieziai/mafs_HVgp_nodup.png", width=15, height=6)

calls_out %>%
  ggplot(aes(x=w, y=medDP)) +
  geom_line(aes(col=pool), alpha=0.7) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggsave("~/Documents/mieziai/dps_HVgp_nodup.png", width=15, height=6)

# Same over genetic map positions
ggplot(calls_out_cm, aes(x=w_cM, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out_cm, aes(x=w_cM, y=medDP)) +
  geom_line(aes(col=pool), alpha=0.7) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

ggplot(calls_out_cm, aes(x=w_cM, y=meanMAF, col=pool)) +
  geom_line(alpha=0.3) +
  geom_smooth(se=F) + 
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggsave("~/Documents/mieziai/mafs_HVgp_nodup_cM.png", width=15, height=6)


ggplot(calls_out, aes(x=w, y=mean11)) +
  geom_line(aes(col=pool)) +
  # scale_y_log10()+
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out, aes(x=w, y=sum11)) +
  geom_line(aes(col=pool)) +
  scale_y_log10()+
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

ggplot(calls_out_cm, aes(x=w_cM, y=sum11+1, col=pool)) +
  geom_line(alpha=0.7) +
  # scale_y_log10() +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  coord_cartesian(ylim=c(0,30000)) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out_cm, aes(x=w_cM, y=mean11, col=pool)) +
  geom_line(alpha=0.3) +
  geom_smooth(se=F) + 
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggsave("~/Documents/mieziai/mean11_HVgp_nodup_cM.png", width=15, height=6)


ggplot(calls_out, aes(x=w, y=meanIndex)) +
  geom_line(alpha=0.3) +
  geom_smooth(se=F) + 
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out_cm, aes(x=w_cM, y=meanIndex)) +
  geom_line(alpha=0.3) +
  geom_smooth(se=F) + 
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")


# ------------------
# chr 4 tests
chr = 1
calls = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                  "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                  .id="pool")
calls = calls[,c("pool", "POS", "REF", "ALT", "DP", "MAF")]

calls = mutate(calls, w = POS %/% 1e6)

dp_calls = filter(calls, pool=="AII")
dp_calls = mutate(dp_calls, kb=POS %/% 1e5 * 100) %>%
  group_by(kb) %>%
  summarize(DP=mean(DP), n=n(), meanMAF=mean(MAF))


# raw depths
dpA = read.table("../depth/dpkb_AII_HVgp_nodups_chr3.txt")
colnames(dpA) = c("kb", "DP")
dpA = mutate(dpA, kb = kb %/% 100 * 100) %>% # 0.1 Mbp windows
  group_by(kb) %>%
  summarize(DP=sum(DP))

p1 = bind_rows("full"=mutate(dpA, DP=DP/1e5), "calls"=dp_calls, .id="source") %>%
  ggplot(aes(x=kb, y=DP+1, col=source)) +
  geom_line(alpha=0.6) + theme_bw() + xlab(NULL) +
  scale_y_log10() + theme(legend.position="top", axis.ticks.x = element_blank(),
                          axis.text.x = element_blank()) 
p2 = ggplot(dp_calls, aes(x=kb, y=n)) +
  geom_line(alpha=0.6, col="grey60") + theme_bw() +
  ylab("n calls/kb") +xlab(NULL) +
  scale_y_log10() + theme(axis.ticks.x = element_blank(),
                          axis.text.x = element_blank())
p3 = ggplot(dp_calls, aes(x=kb, y=meanMAF)) +
  geom_line(alpha=0.6, col="purple") + theme_bw() +
  ylab("mean altAF in calls")
cowplot::plot_grid(p1, p2, p3, nrow=3, rel_heights = c(3,2,2))


bind_rows("full"=mutate(dpA, DP=DP/1e5), "calls"=dp_calls, .id="source") %>%
  group_by(source, Mb=kb %/% 1000) %>%
  filter(Mb<573) %>%
  summarize(DP=mean(DP)) %>%
  ggplot(aes(x=Mb, y=DP+1, col=source)) +
  geom_line(alpha=0.6) + theme_bw() +
  scale_y_log10()

calls4 = calls
mutate(calls4, w=POS %/% 20e6) %>%
  ggplot(aes(x=MAF, group=w, col=w)) + geom_density() +
  theme_bw() + facet_grid(pool~., scales="free_y")

# TODO this for all chrs
mutate(calls, w=POS %/% 20e6) %>%
  ggplot(aes(x=MAF, group=w, col=w)) + geom_density() +
  coord_cartesian(ylim=c(0,10)) +
  theme_bw() + facet_grid(pool~., scales="free_y")

## look at the weird bits for chr4
calls$chunk = cut(calls$POS, c(70e6, 83e6,  300e6, 320e6,  420e6, 440e6,  480e6, 500e6),
                    labels=c("070h", NA,  "mid", NA,   "420m", NA,  "480l"))
ggplot(calls, aes(x=MAF, col=chunk)) + geom_density() +
  theme_bw() + facet_grid(pool~., scales="free_y")

p1 = ggplot(dp_calls, aes(x=kb, y=meanMAF)) +
  geom_line(alpha=0.6, col="purple") +
  theme_bw() + ylab("mean altAF in calls")
p2 = ggplot(dp_calls, aes(x=cM, y=meanMAF)) +
  geom_line(alpha=0.6, col="purple") +
  theme_bw() + ylab("mean altAF in calls")
plot_grid(p1, p2, nrow=2)

# ----------------------------

chr = 4
callsA = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                  "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                  .id="pool")
callsE = bind_rows("AII"=fread(paste0("exons_AII_HVgp_nodups_chr",chr,".txt")),
                  "tw"=fread(paste0("exons_tw_HVgp_nodups_chr",chr,".txt")),
                  .id="pool")
calls = bind_rows("all"=callsA, "exon"=callsE, .id="type")

calls = calls[,c("pool", "type", "POS", "REF", "ALT", "DP", "MAF")]
calls = mutate(calls, w = POS %/% 3e6 * 3)
calls_sum = calls %>%
  group_by(w_cM = w_cM%/%2*2, pool, type) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            chr=chr, meanMAF=mean(MAF))
ggplot(calls_sum, aes(x=w_cM, y=meanMAF, col=pool)) +
  geom_line() +
  # scale_y_log10() +
  facet_grid(type~.) + theme_bw()

# ------------ SNPS ----------------

snps = fread("../refs/HV_M3_variation.vcf.gz")
snps = mutate(snps, CHR=as.numeric(substr(`#CHROM`, 1,1)))
snps = snps[!is.na(snps$CHR), c("CHR", "POS", "REF", "ALT")]

calls = bind_rows("AII"=fread(paste0("callsf_AII_M3_chr",chr,".txt")),
                  "tw"=fread(paste0("callsf_tw_M3_chr",chr,".txt")),
                  .id="pool")
snps_chr = filter(snps, CHR==chr)
snps_chr$known = T
calls = left_join(calls, snps_chr, by=c("POS", "REF", "ALT"))

calls = mutate(calls, w = POS %/% 1e6)

calls_sum = calls %>%
  group_by(w, pool, known=!is.na(known)) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            chr=chr, meanMAF=mean(MAF))
ggplot(calls_sum, aes(x=w, col=pool, y=meanMAF)) +
  geom_line() +
  facet_grid(known~., scales="free_x") +
  theme_bw() + theme(legend.position="bottom")


callsG = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                  "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                  .id="pool")
callsM = bind_rows("AII"=fread(paste0("callsf_AII_M3_chr",chr,".txt")),
                  "tw"=fread(paste0("callsf_tw_M3_chr",chr,".txt")),
                  .id="pool")
calls = bind_rows("GP"=callsG,"M3"=callsM,.id="refgen")
rm(callsG, callsM, callsA, callsE)

calls = mutate(calls, w = POS %/% 1e6)
calls_sum = calls %>%
  group_by(w, refgen, pool) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            chr=chr, meanMAF=mean(MAF))
ggplot(calls_sum, aes(x=w, y=n, col=pool)) +
  geom_line() +
  facet_grid(refgen~., scales="free_x") +
  theme_bw() + theme(legend.position="bottom")

exons = read.table("../refs/Hordeum_vulgare_goldenpromise.GPv1.54.gff3", comment.char="#",sep="\t", quote="")
chr1end = filter(exons, V1=="contig1", V4>460e6, grepl("description", V9))
chr1end = chr1end$V9
chr1end = sapply(strsplit(chr1end, ";"), "[[", 3)
chr1end = sapply(strsplit(chr1end, "description="), "[[", 2)
exonsM = read.table("../refs/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.54.gff3", comment.char="#",sep="\t", quote="")
ggout = tibble()
for(gg in chr1end){
  ggout = bind_rows(ggout, filter(exonsM, grepl(gg, V9, fixed=T)) )
}
ggout

filter(exons, grepl("malic enzyme", V9))
filter(exonsM, grepl("malic enzyme", V9))

callsM = read.table("../tmp_malic_M3.pileup", comment.char = "")
callsG = read.table("../tmp_malic_HVgp.pileup", comment.char = "")
colnames(callsM) = c("CHR", "POS", "REF", "DP.a", "RD.a", "Q.a", "DP.t", "RD.t", "Q.t")
colnames(callsG) = colnames(callsM)
callsM = callsM[250:nrow(callsM),]
callsM$relPOS = 1:nrow(callsM)
callsG$relPOS = 1:nrow(callsG)
tmpM = filter(calls, refgen=="M3", POS>=min(callsM$POS), POS<=max(callsM$POS))
tmpG = filter(calls, refgen=="GP", POS>=min(callsG$POS), POS<=max(callsG$POS)) # 0 calls
tmpM = mutate(tmpM, relPOS = POS-callsM$POS[1])
tmpG = mutate(tmpG, relPOS = POS-callsG$POS[1])

calls_sum = calls %>%
  filter((refgen=="M3" & POS>=min(callsM$POS)-20e3 & POS<=max(callsM$POS)+20e3) |
           refgen=="GP" & POS>=min(callsG$POS)-20e3 & POS<=max(callsG$POS)+20e3) %>%
  group_by(w=POS %/% 1e3, refgen, pool) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            meanMAF=mean(MAF))
ggplot(calls_sum, aes(x=w, y=n, col=pool)) +
  geom_line() + geom_point() +
  facet_wrap(~refgen, scales="free_x") +
  theme_bw() + theme(legend.position="bottom")

tmpG = filter(calls, refgen=="GP", POS %/% 1e3==292721) %>%
  mutate(relPOS = POS-callsG$POS[1])
tmpM = filter(calls, refgen=="M3") %>%
  mutate(relPOS = POS-callsM$POS[1]) %>%
  filter(relPOS > -13000, relPOS < -12000)

left_join(select(callsG, -starts_with("Q")),
          select(callsM, -starts_with("Q")),
          by="relPOS", suffix=c("GP", "M3")) %>%
  select(-starts_with("CHR")) %>%
  filter(grepl("[[:alpha:]]", RD.aGP) | grepl("[[:alpha:]]", RD.tGP) |
           grepl("[[:alpha:]]", RD.aM3) | grepl("[[:alpha:]]", RD.tM3)) %>%
  # filter(tolower(REFGP)!=tolower(REFM3)) %>%
  View

# --------- barke reference
calls = fread(paste0("callsf_twL4_barke.txt"))
calls = calls[,c("CHR", "POS", "REF", "ALT", "DP", "MAF")]
calls = filter(calls, DP>3)

map = read.table("../refs/WGS_ANC.TXT", h=T)
# Most contigs are missing ://
calls = inner_join(calls, map, by=c("CHR"="sequence_id"))
calls = calls[,c("chr", "Mb", "cM", "DP", "MAF")]

calls = mutate(calls, w = Mb %/% 1e6, w_cM = cM %/% 1)
calls_sum = calls %>%
  group_by(w, chr) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            meanMAF=mean(MAF))
calls_sum_cm = calls %>%
  group_by(w_cM, chr) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            meanMAF=mean(MAF))
ggplot(calls_sum_cm, aes(x=w_cM, y=meanMAF)) +
  geom_line(alpha=0.4) +
  geom_smooth(se=F) + 
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
contigs = read.table("../refs/Barke_contig_stats.txt")
map = inner_join(map, contigs, by=c("sequence_id"="V1"))
sum(map$V3/1e6) # only covers 360 Mbp = 10 % of GP1
group_by(map, chr, w=Mb %/% 1e6) %>%
  summarize(cover=sum(V3/1e6)) %>%
  ggplot(aes(x=w, y=cover)) + facet_wrap(~chr, scales="free_x") + geom_line()
group_by(map, chr, w=cM %/% 1) %>%
  summarize(cover=sum(V3/1e6)) %>%
  ggplot(aes(x=w, y=cover)) + facet_wrap(~chr, scales="free_x") + geom_line()

map = read.table("../refs/WGS_ANC.TXT", h=T)
map$strain = sapply(strsplit(map$sequence_id, "_"), "[[", 1)
mapGP = read.table("../refs/barley_morex_x_barke_high_resolution_gbs_map.tsv", h=T)
mapGP = mapGP[,2:4]
mapGP$chr = as.numeric(substr(mapGP$chr, 4, 4))
colnames(mapGP) = c("chr", "Mb", "cM")
mapGP$strain = "GP"
map = bind_rows(map, mapGP)
map = group_by(map, chr, strain) %>%
  mutate(min_old=min(Mb), len_old=max(Mb)-min(Mb)) %>%
  mutate(Mb_new=(Mb-min_old) / len_old)
filter(map, chr==1) %>%
  ggplot(aes(x=Mb_new/1e6, y=cM, col=strain)) +
  geom_point() + theme_bw()

# --------------------
# snpEff annotations

ann = fread("../annots/annot_tw_mono.txt")
nrow(ann)
# already filtered by GT, MQ, DP

# filter strictly by MAF:
ann = filter(ann, MAF>0.9)
nrow(ann)

ann = filter(ann, grepl("HIGH", ANN) | grepl("MODERATE", ANN) | grepl("LOW", ANN))
nrow(ann)

# read in positions of snps called in controls
contr = tibble()
for(chr in 1:7){
  tmp = fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt"))[,1:3]  
  tmp$chr = chr
  contr = bind_rows(contr, tmp)
}
contr = apply(contr, 1, paste, collapse="_")
ann$key = apply(ann[,c("POS", "REF", "ALT", "CHR")], 1, paste, collapse="_")
ann$incontrs = ann$key %in% contr

ann$conseq = ifelse(grepl("HIGH", ann$ANN), "HI", ifelse(grepl("MODERATE", ann$ANN), "MOD", "lo"))


group_by(ann, CHR, incontrs, w=POS %/% 1e6) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=w, y=n)) + 
  geom_point() +
  facet_grid(incontrs~CHR, scales="free_x") +
  theme_bw() + theme(legend.position="bottom")

group_by(ann, CHR, w=POS %/% 1e6) %>%
  summarize(n=n(), meannew=mean(!incontrs)) %>%
  ggplot(aes(x=w, y=meannew)) + 
  geom_point() +
  facet_wrap(~CHR, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")


# w/ genetic map
map = read.table("../refs/WGS_ANC.TXT", h=T)
map$strain = sapply(strsplit(map$sequence_id, "_"), "[[", 1)
map = filter(map, strain=="morex")
# not perfect chr lengths from high-effect snps, but within ~ 5 Mb
chrlens = group_by(ann, CHR) %>%
  summarize(max_new = max(POS))
map = inner_join(map, chrlens, by=c("chr"="CHR"))
map = group_by(map, chr) %>%
  mutate(len_old=max(Mb)) %>%
  mutate(Mb_new=Mb / len_old*max_new)

ann$cM = NA
for(chr in 1:7){
  map_chr = map[map$chr==chr,]
  ann$cM[ann$CHR==chr] = approx(map_chr$Mb_new, map_chr$cM, ann$POS[ann$CHR==chr], rule=2, ties=mean)$y  
}

group_by(ann, CHR, incontrs, w_cM=cM %/% 1) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=w_cM, y=n)) + 
  geom_point() +
  facet_grid(incontrs~CHR, scales="free_x") +
  theme_bw() + theme(legend.position="bottom")

group_by(ann, CHR, w_cM=cM %/% 1) %>%
  summarize(n=n(), meannew=mean(!incontrs)) %>%
  ggplot(aes(x=w_cM, y=meannew)) + 
  geom_point() +
  facet_wrap(~CHR, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

filter(ann, !incontrs) %>%
  mutate(conseq=factor(conseq, levels=c("lo", "MOD", "HI"))) %>%
  ggplot(aes(x=POS/1e6, y=conseq)) + 
  geom_point() +
  facet_wrap(~CHR, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

filter(ann, !incontrs, conseq=="HI") %>% View

# see if we can connect these to any genes w/ names
descr = read.table("../refs/Hordeum_vulgare_goldenpromise.GPv1.54.gff3", comment.char="#",sep="\t", quote="")
descr = filter(descr, grepl("description", V9))
descr = descr$V9
descr = strsplit(descr, ";", fixed=T)
descr_df = data.frame(gene_id = sub("gene_id=", "", sapply(descr, "[[", 4)),
                      description = sub("description=", "", sapply(descr, "[[", 3)))

bads = filter(ann, !incontrs, conseq=="HI")
bads = separate_rows(bads, ANN, sep=",") %>%
  filter(grepl("HIGH", ANN) | grepl("MODERATE", ANN) | grepl("LOW", ANN))

head(bads$ANN)
bads_subfields = strsplit(bads$ANN, "|", fixed=T)
table(sapply(bads_subfields, "[[", 4)==sapply(bads_subfields, "[[", 5))
bads$GENE = sapply(bads_subfields, "[[", 4)
View(bads)
bads = left_join(bads, descr_df, by=c("GENE"="gene_id"))
# nothing special really, ~0.5 % have descrs, and here 1/200

filter(bads, DP4>10, MAF>0.95) %>% View
pairs = paste(bads$REF, bads$ALT, sep="_")
pairs = data.frame(pairs)
# match patterns like CA_CAA, CA_CAAA
pairs$ins = grepl("^([ACGT])([ACGT]+)_\\1\\2([ACGT])(\\3)?$", pairs$pairs)
pairs$del = grepl("^([ACGT])([ACGT]+)([ACGT])(\\3)?_\\1\\2$", pairs$pairs)
pairs$parsed = ifelse(pairs$ins, sub("([ACGT])([ACGT]+)_\\1\\2([ACGT])(\\3)?", "root:\\1\\2 mut:\\3", pairs$pairs),
                ifelse(pairs$del, sub("([ACGT])([ACGT]+)([ACGT])(\\3)?_\\1\\2", "root:\\1\\2 mut:\\3", pairs$pairs), ""))
pairs = separate(pairs, parsed, c("root", "mut"), sep=" ")
pairs$root = sub("root:", "", pairs$root)
pairs$mut = sub("mut:", "", pairs$mut)
pairs$repeated = substring(pairs$root, nchar(pairs$root))==pairs$mut | 
  strrep(substring(pairs$root, nchar(pairs$root)), 2)==pairs$mut
tail(pairs)
bads$repeated = pairs$repeated
# all chrs have some repeats, bit more in chr5 & chr6
filter(bads, is.na(repeated), !duplicated(paste(CHR, POS, ALT))) %>% View


# ----------------
# OLD / NON-POOLED

# load genetic map
map = read.table("../refs/barley_morex_x_barke_high_resolution_gbs_map.tsv", h=T)
map = map[,2:4]
map$chr = as.numeric(substr(map$chr, 4, 4))

calls_out = tibble()
calls_out_cm = tibble()
for(chr in 1:7){
  cold = bind_rows("AII"=fread(paste0("../old/newcalls/callsf_AII-11_HVgp_chr",chr,".txt")),
                   "tw"=fread(paste0("../old/newcalls/callsf_tw-12_HVgp_chr",chr,".txt")),
                   .id="pool")
  cold = mutate(cold, w = POS %/% 2e6 * 2)  # in 2 Mbp windows
  
  # merge in the genetic coords
  map_chr = map[map$chr==chr,]
  range_old = range(map_chr$pos)
  range_new = range(cold$POS)
  map_chr$pos = (map_chr$pos-range_old[1]) / diff(range_old) * diff(range_new) + range_new[1]
  cold$cM = approx(map_chr$pos, map_chr$gbs_cM, cold$POS, rule=2)$y
  cold = mutate(cold, w_cM = cM %/% 1 * 1)  # in 1 cM windows
  
  calls_sum = cold %>%
    group_by(w, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))
  calls_sum_cm = cold %>%
    group_by(w_cM, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))
  calls_out = bind_rows(calls_out, calls_sum)
  calls_out_cm = bind_rows(calls_out_cm, calls_sum_cm)
}

ggplot(calls_out, aes(x=w, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out_cm, aes(x=w_cM, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

chr=5
cold = bind_rows("AII"=fread(paste0("../old/newcalls/callsf_AII-11_HVgp_chr",chr,".txt")),
                 "tw"=fread(paste0("../old/newcalls/callsf_tw-12_HVgp_chr",chr,".txt")),
                 .id="pool")
calls = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                 "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                 .id="pool")
calls_sum = bind_rows("old"=cold, "pool"=calls, .id="proj") %>% 
  mutate(w = POS %/% 2e6 * 2) %>%
  group_by(proj, w, pool) %>%
  summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
            chr=chr, meanMAF=mean(MAF))
ggplot(calls_sum, aes(x=w, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~proj, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")

tmp = inner_join(cold, calls, by=c("pool", "POS", "REF", "ALT"))
tmp = filter(tmp, POS<80e6, POS>50e6)
ggplot(tmp, aes(x=MAF.x, y=MAF.y, col=pool)) +
  geom_point(size=0.7) + theme_bw() +
  xlab("MAF indiv") + ylab("MAF pool")

filter(tmp, pool=="tw") %>%
  ggplot(aes(x=MAF.x, y=MAF.y, col=log(DP.x))) +
  geom_point(size=0.7) + theme_bw() +
  xlab("MAF indiv") + ylab("MAF pool")

filter(tmp, pool=="tw") %>%
  ggplot(aes(x=DP.x, y=DP.y)) +
  scale_y_log10() + scale_x_log10() +
  stat_smooth(method="lm", col="grey70") +
  geom_point(size=0.7) + theme_bw() +
  xlab("DP indiv") + ylab("DP pool")
filter(tmp, pool=="tw", DP.x<120, DP.y<150) %>%
  ggplot(aes(x=DP.x, y=DP.y)) +
  stat_smooth(method="lm", col="grey70") +
  geom_vline(xintercept=median(tmp$DP.x), col="grey50") +
  geom_hline(yintercept=median(tmp$DP.y), col="grey50") +
  geom_point(size=0.7) + theme_bw() +
  xlab("DP indiv") + ylab("DP pool")

# MAFs and even DPs are too similar between pooled and individual runs.
# WTFFFFFFFFFFFFFFF
# simulation:
simn = nrow(tmp)*5
simdp = 50
#simfreq = runif(simn)
#simfreq = rbeta(simn, 0.5, 0.5)
simfreq = rbeta(simn, 1, 5)
simgt = rbinom(simn, 2, simfreq)
simindiv = rbinom(simn, simdp, simgt/2)
simpool = rbinom(simn, simdp, simfreq)
simdf = data.frame(simindiv, simpool) %>%
  filter(simindiv>simdp*0.2, simpool>simdp*0.2)
ggplot(simdf, aes(x=simindiv, y=simpool)) +
  geom_jitter(size=0.7) + theme_bw()

# Trying to see if there are some markers unaffected by this:

# Not just deletions:
filter(tmp, pool=="tw", DP.x<120, DP.y<150) %>%
  mutate(indel=nchar(REF)>1 | nchar(ALT)>1) %>%
  ggplot(aes(x=DP.x, y=DP.y, col=indel)) +
  stat_smooth(method="lm", alpha=0.5) +
  geom_vline(xintercept=median(tmp$DP.x), col="grey50") +
  geom_hline(yintercept=median(tmp$DP.y), col="grey50") +
  geom_point(size=0.7) + theme_bw() +
  scale_color_manual(values=c("skyblue", "orange")) + 
  xlab("DP indiv") + ylab("DP pool")

exons = fread(paste0("exons_tw_HVgp_nodups_chr",chr,".txt"))
tmp$inexons = paste(tmp$POS, tmp$REF, tmp$ALT, sep="_") %in%
  paste(exons$POS, exons$REF, exons$ALT, sep="_")
filter(tmp, pool=="tw", DP.x<120, DP.y<150) %>%
  ggplot(aes(x=DP.x, y=DP.y, col=inexons)) +
  stat_smooth(method="lm", alpha=0.5) +
  geom_vline(xintercept=median(tmp$DP.x), col="grey50") +
  geom_hline(yintercept=median(tmp$DP.y), col="grey50") +
  geom_point(size=0.7) + theme_bw() +
  scale_color_manual(values=c("skyblue", "orange")) + 
  xlab("DP indiv") + ylab("DP pool")

repeats = fread(paste0("../old/newcalls/callsf_tw-12_HVgp_repmasked_chr6.txt"))
repeats$isrepeat = grepl("[a-z]", repeats$REF) | grepl("[a-z]", repeats$ALT)
table(repeats$isrepeat)   # not a lot of the marked repeats though
filter(tmp, pool=="tw", DP.x<120, DP.y<150) %>%
  left_join(repeats[,c("POS", "isrepeat")], by="POS") %>%
  ggplot(aes(x=DP.x, y=DP.y, col=isrepeat)) +
  stat_smooth(method="lm", alpha=0.5) +
  geom_vline(xintercept=median(tmp$DP.x), col="grey50") +
  geom_hline(yintercept=median(tmp$DP.y), col="grey50") +
  geom_point(size=0.7) + theme_bw() +
  scale_color_manual(values=c("skyblue", "orange")) + 
  xlab("DP indiv") + ylab("DP pool")

# will RAW DEPTH be better??
dp.m = fread("../depth/dp_tw-12_HVgp_chr6_mpileup.txt")
dp.m.nb = fread("../depth/dp_tw-12_HVgp_chr6_mpileup_nobaq.txt")
dp.pool = fread("../depth/dp_tw_HVgp_nodups_chr6.txt")
dp.ind = fread("../depth/dp_tw-12_HVgp_chr6.txt")
dp = data.frame("pos"=1:nrow(dp.pool), "pooled"=dp.pool$V1, "indiv"=dp.ind$V1)
colnames(dp.m) = c("pos", "mpileup")
dp.m$mpileupNB = dp.m.nb$V2
dp = left_join(dp, dp.m, by="pos")

# BAQ calculation doesn't matter
filter(dp, pos>18500e3, pos<18520e3) %>%
  ggplot(aes(x=mpileup, y=mpileupNB)) +
  geom_jitter(size=0.5) + 
  theme_bw()

# mpileup vs depth differs slightly, not much
filter(dp, pos>18500e3, pos<18520e3) %>%
  ggplot(aes(x=mpileup, y=indiv)) +
  geom_jitter(size=0.5) + 
  theme_bw()

# differences in depth vs mpileup are mostly due to -A flag in CNVs it seems 
# (not properly paired), which again points to SVs

dp_sum = gather(dp, key="type", value="dp", pooled:mpileup) %>%
  mutate(w=pos %/% 1000) %>%
  group_by(type, w) %>%
  summarize(dp=mean(dp, na.rm=T))
ggplot(dp_sum, aes(x=w, y=dp+1, col=type)) +
  scale_y_log10() +
  geom_point(size=0.5, alpha=0.3) + theme_bw()

dp_sum = filter(dp, pos>18500e3, pos<18800e3) %>%
  gather(key="type", value="dp", pooled:mpileup) %>%
  mutate(w=pos %/% 100 * 0.1) %>%
  group_by(type, w) %>%
  summarize(dp=mean(dp, na.rm=T))
ggplot(dp_sum, aes(x=w, y=dp+1, col=type)) +
  scale_y_log10() + xlab("kbp") +
  coord_cartesian(xlim=c(18500, 18540), ylim=c(5, 90)) +
  geom_line(size=0.5) + theme_bw()

spread(dp_sum, key="type", dp) %>%
  ggplot(aes(x=indiv, y=pooled)) +
  geom_point(size=0.5) + theme_bw()

tmp = inner_join(cold, calls, by=c("pool", "POS", "REF", "ALT"))
tmp = filter(tmp, pool=="tw")

dp$variant = dp$pos %in% tmp$POS
filter(dp, pos>18500e3, pos<18700e3, pos %% 30==0) %>%
  ggplot(aes(x=indiv, y=pooled)) +
  geom_jitter(size=0.5) + theme_bw()

# !! variants seem to concentrate in positions that have
# high DPs in both runs - possibly CNVs?
filter(dp, pos>18500e3, pos<18550e3) %>%
  filter(pooled<120, indiv<150) %>%
  ggplot(aes(x=indiv, y=pooled, col=variant)) +
  geom_point(size=0.7) + theme_bw() +
  scale_color_manual(values=c("skyblue", "orange")) + 
  xlab("DP indiv") + ylab("DP pool")
filter(dp, pos>16000e3, pos<16300e3, variant | pos %%30==0) %>%
  filter(pooled<120, indiv<150) %>%
  ggplot(aes(x=indiv, y=pooled, col=variant)) +
  geom_point(size=0.7) + theme_bw() +
  scale_color_manual(values=c("grey80", "coral")) + 
  xlab("DP indiv") + ylab("DP pool")

filter(dp, pos>14000e3, pos<15000e3, variant | pos %%100==0) %>%
  filter(pooled<120, indiv<150) %>%
  arrange(variant) %>%
  ggplot(aes(x=indiv, y=pooled, col=variant, size=variant)) +
  geom_smooth(method="lm") +
  geom_jitter(width=0.3, height=0.3, alpha=0.6) + theme_bw() +
  scale_size_manual(values=c(0.7, 1.5)) +
  scale_color_manual(values=c("grey80", "coral")) + 
  xlab("DP indiv") + ylab("DP pool")

# more neatly, counts & frac of variable positions over large pieces:
median_pool = median(dp$pooled)
median_indiv = median(dp$indiv)
dp_sum = filter(dp, pos>0, pos<20e6) %>%
  filter(pooled<150, indiv<120, pooled>0, indiv>0) %>%
  group_by(pooled=pooled %/% 3*3, indiv = indiv %/% 3*3) %>%
  summarize(f_var = mean(variant), n=n())
filter(dp_sum, n>50) %>%
  ggplot(aes(x=indiv, y=pooled, col=f_var, size=n)) +
  geom_vline(xintercept=median_indiv, col="grey70") +
  geom_hline(yintercept=median_pool, col="grey70") +
  geom_point() +
  scale_color_distiller(palette="YlGnBu", direction=1) +
  theme_bw() + theme(panel.grid = element_blank()) +
  xlab("DP indiv") + ylab("DP pool")
mutate(dp_sum, n_var = f_var*n) %>%
  ggplot(aes(x=indiv, y=pooled, fill=n_var)) +
  geom_vline(xintercept=median_indiv, col="grey70") +
  geom_hline(yintercept=median_pool, col="grey70") +
  geom_raster(alpha=0.8) +
  scale_fill_distiller(palette="YlGnBu", direction=1) +
  theme_bw() + theme(panel.grid = element_blank()) +
  xlab("DP indiv") + ylab("DP pool")



# TODO:
# A. run SV detector. maybe it is possible to get locations
#  of CNVs and subsequently delete them.
# B. apply a strict "box" boundary to use only the variants near median DP.

# B: corr much weaker if we do this:
filter(dp, pos>15000e3, pos<16000e3, variant | pos %%100==0) %>%
  filter(pooled<120, indiv<150) %>%
  arrange(variant) %>%
  ggplot(aes(x=indiv, y=pooled, col=variant, size=variant)) +
  geom_smooth(method="lm") +
  geom_jitter(width=0.3, height=0.3, alpha=0.6) + theme_bw() +
  geom_vline(data=tibble(), aes(xintercept=c(12, 32)), col="grey50") +
  geom_hline(data=tibble(), aes(yintercept=c(17, 48)), col="grey50") +
  scale_size_manual(values=c(0.7, 1.5)) +
  scale_color_manual(values=c("grey80", "coral")) + 
  xlab("DP indiv") + ylab("DP pool")

filter(dp, pos>15000e3, pos<16000e3, variant | pos %%100==0) %>%
  filter(indiv>12, indiv<32, pooled>17, pooled<48) %>%
  arrange(variant) %>%
  ggplot(aes(x=indiv, y=pooled, col=variant, size=variant)) +
  geom_smooth(method="lm") +
  geom_jitter(width=0.3, height=0.3, alpha=0.6) + theme_bw() +
  geom_vline(data=tibble(), aes(xintercept=c(12, 32)), col="grey50") +
  geom_hline(data=tibble(), aes(yintercept=c(17, 48)), col="grey50") +
  scale_size_manual(values=c(0.7, 1.5)) +
  scale_color_manual(values=c("grey80", "coral")) + 
  xlab("DP indiv") + ylab("DP pool")


# -------------------------
# try plotting MAF of some selected positions
# e.g. by near 0.5 MAF

calls_out = tibble()
calls_out_cm = tibble()
for(chr in 1:7){
  calls = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                    "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                    .id="pool")
  contr_old = fread(paste0("../old/newcalls/callsf_AII-11_HVgp_chr",chr,".txt"))
  contr_old = filter(contr_old, (MAF<0.6 & MAF>0.4 | MAF<0.1), DP<40, DP>10)
  calls = semi_join(calls, contr_old, by=c("POS", "REF", "ALT"))
  
  calls = mutate(calls, w = POS %/% 2e6 * 2)  # in 2 Mbp windows
  
  # merge in the genetic coords
  map_chr = map[map$chr==chr,]
  range_old = range(map_chr$pos)
  range_new = range(calls$POS)
  map_chr$pos = (map_chr$pos-range_old[1]) / diff(range_old) * diff(range_new) + range_new[1]
  calls$cM = approx(map_chr$pos, map_chr$gbs_cM, calls$POS, rule=2)$y
  calls = mutate(calls, w_cM = cM %/% 1 * 1)  # in 1 cM windows
  
  calls_sum = calls %>%
    group_by(w, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))
  calls_sum_cm = calls %>%
    group_by(w_cM, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.9), sum11=sum(MAF>0.9), medDP=median(DP),
              chr=chr, meanMAF=mean(MAF))
  calls_out = bind_rows(calls_out, calls_sum)
  calls_out_cm = bind_rows(calls_out_cm, calls_sum_cm)
}

ggplot(calls_out, aes(x=w, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")
ggplot(calls_out_cm, aes(x=w_cM, y=meanMAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x", nrow=2) +
  theme_bw() + theme(legend.position="bottom")


# ---------------
# SVs from the individual pools
svsT = read.table(gzfile("../old/newcalls/svs/tw/diploidSV.vcf.gz"))
colnames(svsT) = c("CHR", "POS", "ID", "REF", "ALT", "QUALT", "FILTER", "ANN", "FMT", "GT")
svsT = separate(svsT, GT, c("GT", "FT", "GQ", "PL", "PR", "SR"), ":") %>%
  separate_rows(ANN, sep=";") %>%
  extract(ANN, c("ANN_NAME", "ANN_VALUE"), regex="(.*)=(.*)")
svsT = filter(svsT, !is.na(ANN_NAME)) %>%
  spread(ANN_NAME, ANN_VALUE)
svsT$CHR = as.numeric(substr(svsT$CHR, 7, 8))

svsA = read.table(gzfile("../old/newcalls/svs/AII/diploidSV.vcf.gz"))
colnames(svsA) = c("CHR", "POS", "ID", "REF", "ALT", "QUALT", "FILTER", "ANN", "FMT", "GT")
svsA = separate(svsA, GT, c("GT", "FT", "GQ", "PL", "PR", "SR"), ":") %>%
  separate_rows(ANN, sep=";") %>%
  extract(ANN, c("ANN_NAME", "ANN_VALUE"), regex="(.*)=(.*)")
svsA = filter(svsA, !is.na(ANN_NAME)) %>%
  spread(ANN_NAME, ANN_VALUE)
svsA$CHR = as.numeric(substr(svsA$CHR, 7, 8))

table(svsT$GT)
table(svsA$GT) # about half are homozygous
table(svsA$POS %in% svsT$POS)
table(svsT$POS %in% svsA$POS)  # most don't overlap, but this is normal when calling separately

svsT = filter(svsT, GT=="1/1")
svsA = filter(svsA, GT=="1/1")

# manta's quality filters:
svsT = filter(svsT, FILTER=="PASS")
svsA = filter(svsA, FILTER=="PASS")
table(svsT$POS %in% svsA$POS)  #much more overlap, so good to filter

table(svsT$SVTYPE)

chrlens = bind_rows(svsA, svsT) %>%
  group_by(CHR) %>%
  summarize(chrstart=max(POS)/1e6) %>%
  arrange() %>%
  mutate(chrstart = cumsum(lag(chrstart, default=0)))

svsT = inner_join(svsT, chrlens, by="CHR") %>%
  mutate(genpos = POS/1e6 + chrstart)
svsA = inner_join(svsA, chrlens, by="CHR") %>%
  mutate(genpos = POS/1e6 + chrstart)

# SV counts:
bind_rows("AII"=svsA, "tw"=svsT, .id="pool") %>%
  group_by(pool, SVTYPE, w = genpos %/% 2 * 2) %>%
  summarize(n=n()) %>%
  ggplot(aes(x=w, y=n, col=pool)) +
  geom_line() +
  scale_x_continuous(breaks=chrlens$chrstart,
                     minor_breaks=scales::pretty_breaks(100),
                     labels=chrlens$CHR, expand=c(0,0.5)) + 
  facet_grid(SVTYPE~., scales="free_y") +
  theme_minimal() + theme(panel.grid.major.x=element_line(colour="grey50"),
                          panel.grid.minor.y=element_blank())

svsT = select(svsT, -one_of(c("FILTER", "FMT", "GT", "FT", "CIPOS", "CIEND", "GQ", "PL")))
svsA = select(svsA, -one_of(c("FILTER", "FMT", "GT", "FT", "CIPOS", "CIEND", "GQ", "PL")))
svsT$SVLEN = as.numeric(svsT$SVLEN)
svsA$SVLEN = as.numeric(svsA$SVLEN)

# one deletion that is definitely not true (100 Mb in chr 2)
svsT = filter(svsT, POS!=210512827)
# remaining lengths are <10 Mbp:
ggplot(svsT, aes(x=abs(SVLEN), fill=SVTYPE)) +
  geom_density(alpha=0.5) +
  theme_bw()
sum(abs(svsT$SVLEN)/1e6, na.rm=T) # covering ~1 % genome (of the known length ones)

# the interesting region chr4:480-500 Mbp is an SV hotspot, even 
# in this control sample which has a high MAF there too
bind_rows("AII"=svsA, "tw"=svsT, .id="pool") %>%
  filter(SVTYPE!="BND", CHR==4) %>%
  ggplot(aes(x=POS/1e6, y=SVTYPE, col=pool)) +
  geom_jitter(size=0.5, pch=1) +
  theme_minimal()
# for other chrs the indel counts match the MAF peaks too.

bndsT = filter(svsT, SVTYPE=="BND") %>%
  select(-one_of(c("END", "SVTYPE", "SVLEN")))
bndsT$ALT = gsub(".*contig([0-9]):([0-9]*).*", "\\1:\\2", bndsT$ALT)
bndsT = separate(bndsT, ALT, c("END_CHR", "END_POS"), ":", convert=T)
bndsT = inner_join(bndsT, rename(chrlens, endchrstart=chrstart), by=c("END_CHR"="CHR")) %>%
  mutate(genendpos = endchrstart+END_POS/1e6)
bndsT$BND_DEPTH = as.numeric(bndsT$BND_DEPTH)
bndsT$MATE_BND_DEPTH = as.numeric(bndsT$MATE_BND_DEPTH)
bndsT = filter(bndsT, BND_DEPTH>5, MATE_BND_DEPTH>5)

bndsT_chr8 = filter(bndsT, CHR==8 | END_CHR==8, CHR!=END_CHR)
bndsT_chr8$genpos = ifelse(bndsT_chr8$CHR!=8, bndsT_chr8$genpos, bndsT_chr8$genendpos)
bndsT_chr8$genendpos = bndsT_chr8$genpos


bndsT %>%
  filter(CHR!=8, END_CHR!=8, CHR!=END_CHR) %>%
  ggplot(aes(x=genpos, xend=genendpos, y=0)) +
  geom_segment(aes(yend=1), alpha=0.3, col="blue") +
  geom_curve(data=filter(bndsT, CHR!=8, END_CHR==CHR),
             aes(yend=0), col="purple") +
  geom_point(data=filter(bndsT, CHR!=8, END_CHR==CHR),
             aes(y=-0.1), col="purple", pch="|") +
  geom_point(data=bndsT_chr8,
             aes(y=-0.2, x=genpos), col="red", pch="|") +
  scale_x_continuous(breaks=chrlens$chrstart,
                     minor_breaks=scales::pretty_breaks(100),
                     labels=chrlens$CHR, expand=c(0,0.5,0,0)) + 
  scale_y_continuous(breaks=c(0,1), labels=c("pre", "post"), name=NULL) +
  theme_minimal() + theme(panel.grid.major.x=element_line(colour="grey50"),
                          panel.grid.minor.y=element_blank())

bndsT %>%
  filter(CHR==4, CHR!=END_CHR, END_CHR!=8) %>%
  ggplot(aes(x=genpos, xend=genendpos, y=0)) +
  geom_segment(aes(yend=1), alpha=0.3, col="blue") +
  geom_curve(data=filter(bndsT, CHR!=8, END_CHR==CHR),
             aes(yend=0), col="purple") +
  geom_point(data=filter(bndsT, CHR!=8, END_CHR==CHR),
             aes(y=-0.1), col="purple", pch="|") +
  geom_point(data=bndsT_chr8,
             aes(y=-0.2, x=genpos), col="red", pch="|") +
  scale_x_continuous(breaks=chrlens$chrstart,
                     minor_breaks=scales::breaks_width(100),
                     labels=chrlens$CHR, expand=c(0,0.5,0,0)) + 
  scale_y_continuous(breaks=c(0,1), labels=c("pre", "post"), name=NULL) +
  theme_minimal() + theme(panel.grid.major.x=element_line(colour="grey50"),
                          panel.grid.minor.y=element_blank())

# ----------
# Calls on reassembled consensus seqs
chr = 5
cons = bind_rows("AII"=fread(paste0("callsconsf_AII_HVgp_chr",chr,".txt")),
                  "tw"=fread(paste0("callsconsf_tw_HVgp_chr",chr,".txt")),
                 .id="pool")
calls = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                   "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                   .id="pool")
calls = filter(calls, DP>5, MQ>30)
cons = filter(cons, DP>5, MQ>30)

calls_conspool = fread(paste0("callsconspoolf_tw_HVgp_chr",chr,".txt"))
calls_conspool$pool = "tw"
calls_conspool = filter(calls_conspool, DP>5, MQ>30)

calls_sum = bind_rows("GP"=calls, "cons"=cons, "cons.pool"=calls_conspool, .id="ref") %>%
  filter(DP>10, DP<60) %>% # this can be included or not, works both ways
  group_by(w=POS %/% 2e6 * 2, pool, ref) %>%
  summarize(meanMAF=mean(MAF), mean11=mean(MAF>0.9), n=n())

ggplot(calls_sum, aes(x=w, y=meanMAF, col=pool)) +
  geom_line() + facet_grid(ref~.) +
  theme_bw()

ggplot(calls_sum, aes(x=w, y=meanMAF, col=pool)) +
  geom_line() + facet_grid(ref~.) +
  annotate(geom="point", y=1, x=c(506.1, 506.9), col="purple", alpha=0.5) +
  theme_bw()

map = read.table("../refs/barley_morex_x_barke_high_resolution_gbs_map.tsv", h=T)
# merge in the genetic coords
map_chr = map[map$chr==paste0("chr",chr,"H"),]
# basic rescaling, to at least match the length of GP chroms:
range_old = range(map_chr$pos)
range_new = range(calls$POS)
map_chr$pos = (map_chr$pos-range_old[1]) / diff(range_old) * diff(range_new) + range_new[1]
calls$cM = approx(map_chr$pos, map_chr$gbs_cM, calls$POS, rule=2)$y
cons$cM = approx(map_chr$pos, map_chr$gbs_cM, cons$POS, rule=2)$y

calls_sum_cm = bind_rows("GP"=calls, "cons"=cons, .id="ref") %>%
  mutate(w_cM = cM %/% 1 * 1) %>%
  group_by(w_cM, pool, ref) %>%
  summarize(meanMAF=mean(MAF), mean11=mean(MAF>0.9), n=n())
ggplot(calls_sum_cm, aes(x=w_cM, y=meanMAF, col=pool)) +
  geom_line() + facet_grid(ref~.) +
  theme_bw()

aii = fread("callsf_AII_HVgp_nodups_chr5.txt")
aiiindiv = fread("../old/newcalls/callsf_AII-11_HVgp_chr5.txt")
aii$ID = paste(aii$POS, aii$REF, aii$ALT, sep="_")
aiiindiv$ID = paste(aiiindiv$POS, aiiindiv$REF, aiiindiv$ALT, sep="_")

table(aii$ID %in% aiiindiv$ID)
table(aiiindiv$ID %in% aii$ID[aii$GT=="1/1"])
table(aii$ID[aii$GT=="1/1"] %in% aiiindiv$ID[aiiindiv$GT=="1/1"])
