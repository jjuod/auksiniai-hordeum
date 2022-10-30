library(dplyr)
setwd("~/Documents/mieziai/calls/")

# prepare gene lists
refM3 = read.table("../refs/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.54.gff3", quote="", h=F, sep="\t")
refM3 = filter(refM3, V3=="gene")

refM3 = mutate(refM3, id=paste(V1, paste(V4, V5, sep="-"), sep=":"))
refM3$gene.M3 = gsub(".*gene_id=(.*);.*", "\\1", refM3$V9)

refGP = read.table("../refs/Hordeum_vulgare_goldenpromise.GPv1.54.gff3", quote="", h=F, sep="\t")
refGP = filter(refGP, V3=="gene")

refGP = mutate(refGP, id=paste(V1, paste(V4, V5, sep="-"), sep=":"))
refGP$gene.GP = gsub(".*gene_id=(.*);.*", "\\1", refGP$V9)

# ---------------------------
# attach both gene names to results

res_hi = read.table("genes/blastres_regions_hi.txt")
colnames(res_hi) = c("qregion", "qlen", "qstart", "qend", "schrom", "sstart", "send", "eval", "slen", "pident")
res_mod = read.table("genes/blastres_regions_mod.txt")
colnames(res_mod) = c("qregion", "qlen", "qstart", "qend", "schrom", "sstart", "send", "eval", "slen", "pident")

# select best match by eval, and when it is tied (at 0), take by %identity
res_hi = group_by(res_hi, qregion) %>%
  top_n(1, -eval) %>%
  filter(slen>0.3*qlen) %>%
  top_n(1, pident)
res_hi = left_join(res_hi, refGP[,c("id", "gene.GP")], by=c("qregion"="id"))

res_mod = group_by(res_mod, qregion) %>%
  top_n(1, -eval) %>%
  filter(slen>0.2*qlen) %>%
  slice_max(n=1, order_by=pident, with_ties=F)
res_mod = left_join(res_mod, refGP[,c("id", "gene.GP")], by=c("qregion"="id"))

all(res_mod$schrom=="5H")
all(res_hi$schrom=="5H")

find_match_M3 = function(df){
  df$match.M3 = NA
  tmp_M3 = filter(refM3, V1=="5H")
  for(i in 1:nrow(df)){
    # find any overlapping genes
    allmatches = filter(tmp_M3, V4<df$send[i], V5>df$sstart[i])
    if(nrow(allmatches)>0) df$match.M3[i] = paste(allmatches$gene.M3, collapse=",")
  }
  return(df)
}
res_hi2 = find_match_M3(res_hi) %>%
  unite("spos", sstart:send, sep="-") %>%
  unite("match.region", schrom:spos, sep=":")
res_mod2 = find_match_M3(res_mod) %>%
  unite("spos", sstart:send, sep="-") %>%
  unite("match.region", schrom:spos, sep=":")

# attach conversions to annotated variants
hits_hi = read.table("final_tw_chr5_hi.csv", sep="\t", h=T)
hits_hi = left_join(hits_hi, res_hi2, by=c("GeneName"="gene.GP"))
  
hits_mod = read.table("final_tw_chr5_mod.csv", sep="\t", h=T)
hits_mod = left_join(hits_mod, res_mod2, by=c("GeneName"="gene.GP"))

write.table(hits_hi, "finalwblast_tw_chr5_hi.csv", quote=F, sep="\t", row.names=F)
write.table(hits_mod, "finalwblast_tw_chr5_mod.csv", quote=F, sep="\t", row.names=F)
