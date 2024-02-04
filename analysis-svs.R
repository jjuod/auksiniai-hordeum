library(tidyr)
library(dplyr)
library(ggplot2)

setwd("~/Documents/mieziai/2024/svs/")


# control variants, if needed:
# svs_aii = read.csv("dysgu_AII-11_filtered.vcf", sep="\t", comment.char = "#", h=F)

# vcfs already contain only variants unique to the tw-12 individual
# and passing some quality filters
svs_aii = read.csv("dysgu_tw-12_unique_filteredafter.vcf", sep="\t", comment.char = "#", h=F)

colnames(svs_aii) = c("CHR", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GT")
nrow(svs_aii)
table(svs_aii$CHR)

# keep only the main chromosome contigs
svs_aii = filter(svs_aii, CHR %in% paste0(1:7, "H"))
nrow(svs_aii)
table(svs_aii$FILTER)
svs_aii = filter(svs_aii, FILTER=="PASS")
nrow(svs_aii)

# parse some fields
svs_aii$END = as.numeric(sub(".*;END=([^;]*);.*", "\\1", svs_aii$INFO))
svs_aii$CHR2 = sub(".*;CHR2=([^;]*);.*", "\\1", svs_aii$INFO)
svs_aii$SVLEN = as.numeric(sub(".*;SVLEN=([^;]*);.*", "\\1", svs_aii$INFO))
svs_aii$SVTYPE = sub(".*;SVTYPE=([^;]*);.*", "\\1", svs_aii$INFO)
# for translocations:
svs_aii$CHR2_POS = as.numeric(sub(".*;CHR2_POS=([^;]*);.*", "\\1", svs_aii$INFO))
svs_aii$END = ifelse(svs_aii$SVTYPE=="TRA", svs_aii$CHR2_POS, svs_aii$END)

table(svs_aii$CHR)

# ------------------
# --- PLOTTING -----

library(circlize)

# keep only larger and homozygous SVs
df = filter(svs_aii, SVLEN>1000, substr(GT, 1,3)=="1/1") %>%
    mutate(POS = POS/ 1e6, END = END/1e6)

dels = filter(df, SVTYPE=="DEL")
ins = filter(df, SVTYPE=="DUP")

circos.initialize(sectors=df$CHR, x=df$POS)
circos.track(df$CHR, ylim=c(0,1),
             panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] + mm_y(5), 
                             CELL_META$sector.index)
                 circos.axis(labels.cex = 0.6)
             })

circos.trackPoints(dels$CHR, dels$POS, rep(0.8, nrow(dels)), pch=17, col="purple")
circos.trackPoints(ins$CHR, ins$POS, rep(0.4, nrow(ins)), pch=17, col="darkblue")
circos.trackHist(df$CHR, df$POS, bin.size=10, col="grey10")
for(i in seq_along(df$CHR)){
    if(df$SVLEN[i]>1e3 & df$SVTYPE[i]!="DEL" & df$SVTYPE[i]!="DUP"){
        circos.link(df$CHR[i], df$POS[i], df$CHR2[i], df$END[i], col="orangered")   
    }
}

svs_aii_fixed = filter(svs_aii, substr(GT, 1,3)=="1/1") %>%
    mutate(POS=POS/1e6, END=END/1e6) %>%
    filter(SVTYPE!="TRA", SVLEN<200e6)

gff = read.table("../../refs/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.54.gff3", comment.char="#", quote="", sep="\t")
gff = filter(gff, V3=="gene")
gff = mutate(gff, POS=V4/1e6, END=V5/1e6)
gff = select(gff, one_of(c("V1", "POS", "END"))) %>% rename(GENESTART=POS, GENEEND=END)

df4 = tibble()
for(chr in paste0(1:7, "H")){
    t1 = filter(svs_aii_fixed, CHR==chr)
    t2 = filter(gff, V1==chr)
    out = crossing(t1, t2) %>% filter(POS<GENEEND & END>GENESTART)
    df4 = bind_rows(df4, out)
}
nrow(df4)

df2 = filter(svs_aii_fixed, CHR=="5H")
df2 = filter(df2, SVLEN>1e3, POS>530)

gff = filter(gff, V1=="5H", V4>530e6)

# check overlap
df3 = crossing(df2, gff) %>%
    filter(POS<GENEEND & END>GENESTART)

ggplot(df4) +
    geom_segment(aes(x=POS, xend=END, y=seq_along(POS)%%10, yend=10, col=SVTYPE),
                    lwd=1) +
    geom_segment(aes(x=POS, xend=END, y=5, yend=5, col=SVTYPE),
                 lwd=10, col="darkblue", data=gff)

# one more circ
png("../result-compare.png", res = 180, height=1000, width=1000)
dels = filter(df4, SVTYPE=="DEL")
ins = filter(df4, SVTYPE=="DUP" | SVTYPE=="INS")
circos.initialize(sectors=df4$CHR, x=df4$POS)
circos.track(df4$CHR, ylim=c(0,1),
             panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] + mm_y(5), 
                             CELL_META$sector.index)
                 circos.axis(labels.cex = 0.6)
             })

circos.trackPoints(dels$CHR, dels$POS, rep(0.8, nrow(dels)), pch=17, col="purple")
circos.trackPoints(ins$CHR, ins$POS, rep(0.4, nrow(ins)), pch=17, col="darkblue")
circos.trackHist(df4$CHR, df4$POS, bin.size=10, col="grey10")
for(i in seq_along(df$CHR)){
    if(df4$SVLEN[i]>1e3 & df4$SVTYPE[i]!="DEL" & df4$SVTYPE[i]!="DUP"){
        circos.link(df4$CHR[i], df4$POS[i], df4$CHR2[i], df4$END[i], col="orangered")   
    }
}
dev.off()

