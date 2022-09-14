library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
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

# zoom in 
calls_sum = filter(calls2, POS>460e6) %>%
  mutate(w=POS %/% 30e3) %>% 
  group_by(w) %>%
  summarize(twMAC=mean(twMAC), AIIMAC=mean(AIIMAC), n=n())

gather(calls_sum, key="pool", value="MAC", twMAC:AIIMAC) %>%
  ggplot(aes(x=w, y=MAC, col=pool)) + 
  geom_line() +
  theme_bw()

filter(calls2, POS>350e6, twMAC==2, AIIMAC==0) %>%
  ggplot(aes(x=POS, y=DP)) +
  geom_point() +
  theme_bw()
filter(calls2, POS>465e6, twMAC==2, AIIMAC==0) %>%
  ggplot(aes(x=POS/1e3, y=DP)) +
  geom_point() +
  theme_bw()

filter(calls2, POS>468e6) %>%
  mutate(w=POS %/% 10e3) %>%
  group_by(w) %>%
  summarize(twMAC=mean(twMAC), AIIMAC=mean(AIIMAC), n=n()) %>%
  gather(key="pool", value="MAC", twMAC:AIIMAC) %>%
  ggplot(aes(x=w, y=MAC, col=pool)) +
  geom_point(size=0.8) +
  theme_bw()


## --------------
# reading in single-sample depths
callsA = fread("calls_AII_HVgp_nodups_chr1.txt")
callsT = fread("calls_tw_HVgp_nodups_chr1.txt")

cts = group_by(callsA, GT) %>%
  summarize(n=n(), meanMapQ=mean(MQ), meanDP=mean(DP)) %>%
  arrange(desc(n))
View(cts)

cts = group_by(callsT, GT) %>%
  summarize(n=n(), meanMapQ=mean(MQ), meanDP=mean(DP)) %>%
  arrange(desc(n))
View(cts)

# still don't like that DP4 is quite filtered out
calls_out = tibble()
for(chr in 1:7){
  callsA = fread(paste0("calls_AII_HVgp_nodups_chr",chr,".txt"))
  callsT = fread(paste0("calls_tw_HVgp_nodups_chr",chr,".txt"))
  
  callsA = separate(callsA, "DP4", c("refF", "refR", "altF", "altR"), sep=",")
  callsT = separate(callsT, "DP4", c("refF", "refR", "altF", "altR"), sep=",")
  callsA = mutate(callsA, ref=as.numeric(refF)+as.numeric(refR),
                  alt=as.numeric(altF)+as.numeric(altR)) %>%
    select(-one_of(c("refF", "refR", "altF", "altR")))
  callsT = mutate(callsT, ref=as.numeric(refF)+as.numeric(refR),
                  alt=as.numeric(altF)+as.numeric(altR)) %>%
    select(-one_of(c("refF", "refR", "altF", "altR")))
  
  calls_sum = bind_rows("AII"=callsA, "tw"=callsT, .id="pool") %>%
    mutate(w=POS %/% 1e6) %>%
    group_by(w, pool) %>%
    summarize(n=n(), MAF=mean(alt/(alt+ref)), chr=chr)
  calls_out = bind_rows(calls_out, calls_sum)
}

ggplot(calls_out, aes(x=w, y=MAF)) +
  geom_line(aes(col=pool)) +
  facet_wrap(~chr, scales="free_x") +
  theme_bw()

## raw depths on a sample
dp_wdup = read.table("../depth/dp_AII_HVgp.txt")
dp_nodup = read.table("../depth/dp_AII_HVgp_nodups.txt")
dp_m3 = read.table("../depth/dp_AII_M3.txt")
dp_tw_wdup = read.table("../depth/dp_tw_HVgp.txt")
dp_tw_nodup = read.table("../depth/dp_tw_HVgp_nodups.txt")
dp_tw_m3 = read.table("../depth/dp_tw_M3.txt")

dp = data.frame(POS=dp_wdup$V2, DP_WDUP=dp_wdup$V3, DP_NODUP=dp_nodup$V3,
                DP_M3=dp_m3$V3, DP_TW_WDUP=dp_tw_wdup$V3,
                DP_TW_NODUP=dp_tw_nodup$V3, DP_TW_M3=dp_tw_m3$V3)
rm(dp_m3, dp_nodup, dp_tw_m3, dp_tw_nodup, dp_tw_wdup, dp_wdup)
dp2 = gather(dp, key="stat", value="cnt", DP_WDUP:DP_TW_M3) %>%
  mutate(ref=ifelse(stat %in% c("DP_M3", "DP_TW_M3"), "M3", "gp"),
         dups=ifelse(stat %in% c("DP_NODUP", "DP_TW_NODUP"), "removed", "kept"),
         pool=ifelse(stat %in% c("DP_TW_WDUP", "DP_TW_NODUP", "DP_TW_M3"), "tw", "AII"))
dp2_sum = group_by(dp2, w=POS %/% 5e3, dups, pool, ref) %>%
  summarize(DP=sum(cnt))

# check duplicate removal:
filter(dp2_sum, ref=="gp") %>%
  ggplot(aes(x=w, y=DP+1, col=pool, lty=dups)) +
  geom_line() + theme_bw() +
  scale_y_log10()
# dup removal doesn't seem to change much, even for the spikes.

dp2 = filter(dp2, dups=="kept")
dp2_sum = filter(dp2_sum, dups=="kept")
dp2_sum %>%
  ggplot(aes(x=w, y=DP+1, col=pool)) +
  facet_grid(ref~.) +
  geom_line() + theme_bw() +
  scale_y_log10()

# genes

g_m3 = read.table("../refs/M3_chr2.gff3", sep="\t")
g_gp = read.table("../refs/gp_chr2.gff3", sep="\t")
g_m3 = filter(g_m3, V4<2e6)
g_gp = filter(g_gp, V4<2e6)
gs = bind_rows("gp"=g_gp, "M3"=g_m3, .id="ref")
gs = mutate(gs, V4=V4 / 5e3, V5=V5 / 5e3) %>%
  filter(V3=="gene")
regexmatches = regexpr("description=[^;]*", gs$V9)
gs$descr = ""
gs$descr[regexmatches!=-1] = regmatches(gs$V9, regexmatches)
# :(


dp2_sum %>%
  ggplot(aes(x=w, y=DP+1, col=pool)) +
  facet_grid(ref~.) +
  geom_line() + theme_bw() +
  geom_segment(data=gs, aes(x=V4, xend=V5, y=1, yend=1, col=descr), size=5) +
  scale_y_log10() + theme(legend.position="bottom")

# ---------------------

# calls / per-sample depths

callsA = fread("calls_AII_HVgp_nodups_chr2.txt")
callsA = filter(callsA, POS<=2e6)
callsT = fread("calls_tw_HVgp_nodups_chr2.txt")
callsT = filter(callsT, POS<=2e6)
callsA = bind_rows("AII"=callsA, "tw"=callsT, .id="pool")

callsA = separate(callsA, "DP4", c("refF", "refR", "altF", "altR"), sep=",")
callsA = mutate(callsA, ref=as.numeric(refF)+as.numeric(refR),
                alt=as.numeric(altF)+as.numeric(altR)) %>%
  select(-one_of(c("refF", "refR", "altF", "altR")))

callsA_sum = group_by(callsA, pool, w=POS %/% 5e3) %>%
  summarize(ncalls=n(), ncalls_dp10=sum(DP>10 & MQ>30),
            altAF=sum(alt)/(sum(ref)+sum(alt)),
            meanMQ=mean(MQ))

sum1 = filter(dp2_sum, ref=="gp") %>%
  left_join(callsA_sum, by=c("w", "pool"))

p1 = sum1 %>% 
  ggplot(aes(x=w, y=DP+1, col=pool)) +
  geom_line() + scale_y_log10() + theme_bw() +
  theme(legend.position = "top")

p2 = sum1 %>% filter(pool=="tw") %>%
  gather(key="dp", value="ncalls", ncalls:ncalls_dp10) %>%
  ggplot(aes(x=w, y=ncalls, lty=dp)) +
  geom_line() + theme_bw() + 
  geom_line(data=filter(sum1, pool=="AII"), aes(y=5*meanMQ), col="skyblue", lty="solid") +
  theme(legend.position="bottom")

p3 = sum1 %>%
  ggplot(aes(x=w, y=altAF, col=pool)) +
  geom_line() + theme_bw() +
  theme(legend.position = "none")

cowplot::plot_grid(p1, p3, p2, nrow=3)

# informative plots:
# DP in one, DP in other
# ncalls in one, ncalls in other
# maf in each??


# ------------- 
callsA = fread("calls_AII_HVgp_nodups_chr2.txt")
callsA = filter(callsA, POS<=350e6)
callsM = fread("tmp_calls_M3_AII_chr2.txt")
colnames(callsM) = colnames(callsA)
callsM = filter(callsM, POS<=350e6)
callsA = bind_rows("M3"=callsA, "gp"=callsM, .id="refgen")

callsA = separate(callsA, "DP4", c("refF", "refR", "altF", "altR"), sep=",")
callsA = mutate(callsA, ref=as.numeric(refF)+as.numeric(refR),
                alt=as.numeric(altF)+as.numeric(altR)) %>%
  select(-one_of(c("refF", "refR", "altF", "altR")))

callsA_sum = group_by(callsA, refgen, w=POS %/% 200e3) %>%
  summarize(ncalls=n(), ncalls_dp10=sum(DP>10 & MQ>30),
            ncalls_mono=sum(ref==0), ncalls_mono1=sum(ref<=1),
            altAF=sum(alt)/(sum(ref)+sum(alt)),
            meanMQ=mean(MQ))

gather(callsA_sum, key="stat", value="y", c(ncalls,ncalls_dp10)) %>%
  ggplot(aes(x=w, y=y+1, col=refgen, lty=stat)) +
  geom_line() + scale_y_log10() +
  theme_bw()
callsA_sum %>%
  ggplot(aes(x=w, y=meanMQ, col=refgen)) +
  geom_line() + theme_bw()
gather(callsA_sum, key="stat", value="y", c(ncalls_mono,ncalls_mono1)) %>%
  ggplot(aes(x=w, y=y, col=refgen, lty=stat)) +
  geom_line() + theme_bw()
# it would seem that the GP reference is actually less similar to ours
# (at least AII)? M3 has fewer calls & less weird regions.
# Depth & MQ filters don't change anything for either.

gather(callsA_sum, key="stat", value="y", ncalls:meanMQ) %>%
  ggplot(aes(x=w, y=y, col=refgen)) +
  geom_line() + 
  facet_grid(stat~., scales="free_y") + theme_bw()
