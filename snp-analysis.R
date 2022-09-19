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

calls = fread("calls_HVgp_nodups_chr1.txt")
cts = group_by(calls, AII, tw) %>%
  summarize(n=n(), meanMapQ=mean(MQ)) %>%
  arrange(desc(n))
View(cts)

# notably much fewer non-ref GTs in AII vs. tw.
# tw have much more inbreeding I guess?

calls = separate(calls, "DP4", c("refF", "refR", "altF", "altR"), sep=",")

# TODO samtools depth or mpileup+vcf call per file?

# for now, select by GT:
calls2 = separate(calls, "tw", c("twL", "twR"), sep="/")
calls2 = separate(calls2, "AII", c("AIIL", "AIIR"), sep="/")
# remove the weird missing haplos  --- TODO look into them
calls2 = filter(calls2, twR!=".", twL!=".", AIIL!=".", AIIR!=".")
# combine all ALT alleles
calls2 = mutate(calls2, twR=pmin(1, as.numeric(twR)),
                twL=pmin(1, as.numeric(twL)),
                AIIL=pmin(1, as.numeric(AIIL)),
                AIIR=pmin(1, as.numeric(AIIR)))
calls2 = mutate(calls2, twMAC=twL+twR, AIIMAC=AIIL+AIIR)

# window
calls_sum = mutate(calls2, w=POS %/% 5e5) %>% 
  group_by(w) %>%
  summarize(twMAC=mean(twMAC), AIIMAC=mean(AIIMAC), n=n())

# snp frequency per kb basically
ggplot(calls_sum, aes(x=w, y=log10(n))) + 
  geom_line() +
  theme_bw()

# ALT af
gather(calls_sum, key="pool", value="MAC", twMAC:AIIMAC) %>%
  ggplot(aes(x=w, y=MAC, col=pool)) + 
  geom_line() +
  theme_bw()


# -----------------------------------
# --------- CLEAN BIT ---------------

# TODO
# informative plots:
# DP in one, DP in other
# ncalls in one, ncalls in other
# maf in each??

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
  calls = bind_rows("AII"=fread(paste0("callsf_AII_HVgp_nodups_chr",chr,".txt")),
                     "tw"=fread(paste0("callsf_tw_HVgp_nodups_chr",chr,".txt")),
                    .id="pool")
  # calls = bind_rows("AII"=fread(paste0("callsf_AII_M3_chr",chr,".txt")),
  #                   "tw"=fread(paste0("callsf_tw_M3_chr",chr,".txt")),
  #                   .id="pool")
  
  # Filtering for mapping quality:
  calls = filter(calls, MQ>30)
  
  # Filtering by depth:
  calls = filter(calls, DP>=15, DP<90)  # about 10 % loss
  
  calls = calls[,c("pool", "POS", "REF", "ALT", "DP", "MAF")]
  
  # Filtering only those present in both samples:
  # calls = inner_join(filter(calls, pool=="AII"),
  #                    filter(calls, pool=="tw"),
  #                    by=c("POS", "REF", "ALT"),
  #                    suffix=c(".AII", ".tw")) %>%
  #   select(-one_of(c("pool.AII", "pool.tw")))
  # # This is important! otherwise depth etc. can really differ
  # # depending on which markers are included
  # 
  # # Filtering to remove roughly monomorphic markers:
  # calls = filter(calls, MAF.AII<0.8 | MAF.tw<0.8)
  # Keep only unique ones for either pool:
  callsBOTH = inner_join(filter(calls, pool=="AII"), filter(calls, pool=="tw"),
                         by=c("POS", "REF", "ALT"))
  calls = anti_join(calls, callsBOTH, by=c("POS", "REF", "ALT"))
  
  calls = mutate(calls, w = POS %/% 2e6 * 2)  # in 2 Mbp windows
  
  # merge in the genetic coords
  map_chr = map[map$chr==chr,]
  # basic rescaling, to at least match the length of GP chroms:
  # TODO
  range_old = range(map_chr$pos)
  range_new = range(calls$POS)
  map_chr$pos = (map_chr$pos-range_old[1]) / diff(range_old) * diff(range_new) + range_new[1]
  calls$cM = approx(map_chr$pos, map_chr$gbs_cM, calls$POS, rule=2)$y
  calls = mutate(calls, w_cM = cM %/% 1 * 1)  # in 1 cM windows
  
  # calls = pivot_longer(calls, c(DP.AII, DP.tw, MAF.AII, MAF.tw),
  #              names_to=c("stat", "pool"), names_sep="\\.", values_to="val") %>%
  #   spread("stat", "val")

  calls_sum = calls %>%
    group_by(w, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.8), sum11=sum(MAF>0.8), medDP=median(DP), chr=chr)
    # summarize(n=n(), meanMAF=mean(MAF), medDP=median(DP), chr=chr)
  calls_sum_cm = calls %>%
    group_by(w_cM, pool) %>%
    summarize(n=n(), mean11=mean(MAF>0.8), sum11=sum(MAF>0.8), medDP=median(DP), chr=chr)
    # summarize(n=n(), meanMAF=mean(MAF), medDP=median(DP), chr=chr)

  calls_out = bind_rows(calls_out, calls_sum)
  calls_out_cm = bind_rows(calls_out_cm, calls_sum_cm)
}

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
  scale_y_log10() +
  facet_wrap(~chr, scales="free_x", nrow=2) +
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
chr = 3
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


