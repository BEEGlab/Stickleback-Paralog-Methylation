##########################################################################################
#R script for final analysis of methylation across individuals. Must run all previous scripts first

library(plyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gplots)
library(OneR) ##for binning


##########################################################################################
meta<-read.table("SraRunTable_Metzger2018.txt",sep="\t",header=T) ##Gets information about samples
meta<-meta[,c(7,1,4,9,10,11,12)]
meta$replicate<-sapply(strsplit(as.character(meta$replicate), " "), `[`, 3)

# load the data AS LOOP
files <- list.files(path=".", pattern="*SRR.*bismark.cov.gz.meth.prom.gene.dups.csv", full.names=TRUE, recursive=TRUE)
# all genes with their promoter information
filesall <- list.files(path=".", pattern="*SRR.*bismark.cov.gz.meth.prom.gene.full.csv", full.names=TRUE, recursive=FALSE)

########################################################################################## 
########################################################################################## 
########################################################################################## 
# Loop through each individual to merge together
allgenes <- vector("list", length(files))
allg<-read.csv(filesall[1],row.names=1)
allg_info<-allg[,names(allg) %in% c("ensembl_gene_id","paralogue_ensembl_gene_id","last_common_ancestor","Cat","ohno","Gene","avg_pct_id","CNV.x","NumInd.x","NumPop.x","average")]
allg<-allg[,names(allg) %in% c("ensembl_gene_id","last_common_ancestor","Cat","ohno","Gene")]

##Save duplicate gene methylation info
data <- vector("list", length(files)) ##need to initiate?
dupg<-read.csv(files[1],row.names=1)
avgs_all<-(dupg[,names(dupg) %in% c("Cat","last_common_ancestor","HK","SP","ovaries","testes","em","lm","mg","eo","X24","average","p_HK","p_SP","p_ovaries","p_testes","p_em","p_lm","p_mg","p_eo","p_X24","p_average")]) ##expression
dupg_info<-dupg[,names(dupg) %in% c("ensembl_gene_id","paralogue_ensembl_gene_id","last_common_ancestor","Cat","ohno","Gene","avg_pct_id","dN","dS","CNVs","NumInds","NumPops","avg_exp","div_exp")]
dupg<-dupg[,names(dupg) %in% c("ensembl_gene_id","paralogue_ensembl_gene_id","last_common_ancestor","Cat","ohno")]

for (f in 1:length(files)){
  ###Get all genes, not just duplicates
  allgenes[[f]]<-read.csv(filesall[f],row.names=1)
  allgenes[[f]]<-allgenes[[f]][,names(allgenes[[f]]) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation")]
  all<-allgenes[[f]]
  colnames(all)<-c("ensembl_gene_id",paste(meta[f,4],"prom",sep="_"),paste(meta[f,4],"prom_counts",sep="_"),paste(meta[f,4],"gb",sep="_"))
  allg<-merge(allg,all,by="ensembl_gene_id")
  
  ###Get all duplicate pairs
  data[[f]]<-read.csv(files[f],row.names=1)
  data[[f]]<-data[[f]][,names(data[[f]]) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation","avg_prom_meth","div_prom_meth","sum_prom_counts","avg_gb_meth","div_gb_meth")]
  dup<-data[[f]]
  colnames(dup)<-c("ensembl_gene_id",paste(meta[f,4],"promoter_methylation",sep="_"),paste(meta[f,4],"prom_counts",sep="_"),paste(meta[f,4],"gb_methylation",sep="_"),paste(meta[f,4],"p_promoter_methylation",sep="_"),paste(meta[f,4],"p_prom_counts",sep="_"),paste(meta[f,4],"p_gb_methylation",sep="_"),paste(meta[f,4],"avg_prom_meth",sep="_"),paste(meta[f,4],"div_prom_meth",sep="_"),paste(meta[f,4],"sum_prom_counts",sep="_"),paste(meta[f,4],"avg_gb_meth",sep="_"),paste(meta[f,4],"div_gb_meth",sep="_"))
  dupg<-merge(dupg,dup,by="ensembl_gene_id")
  
  ###All genes in different format ...
  rownames(allgenes[[f]])<-allgenes[[f]][,1]
  allgenes[[f]]$sample<-meta[f,4]
  allgenes[[f]]$rep<-meta[f,5]
  allgenes[[f]]$sex<-meta[f,6]
  allgenes[[f]]$treatment<-meta[f,7]
  
  ###Duplicate pairs in different format ...
  rownames(data[[f]])<-data[[f]][,1]
  data[[f]]$sample<-meta[f,4]
  data[[f]]$rep<-meta[f,5]
  data[[f]]$sex<-meta[f,6]
  data[[f]]$treatment<-meta[f,7]
}

write.csv(allg, "Combined.gene.methylation.csv",row.names=F)
write.csv(dupg, "Combined.dup.methylation.csv",row.names=F)

########################################################################################## 
########################################################################################## 
########################################################################################## 

## Reorder levels/categories for plotting:
allg_info$Cat <- factor(allg_info$Cat, levels=c("Young","Medium","Old","Singleton"))
allg_info$ohno <- factor(allg_info$ohno, levels=c("SSD","ohnolog","single"))
## LCA -Few Eupercaria and Perciformes (<15), group with Percomorphaceae
allg_info$last_common_ancestor[allg_info$last_common_ancestor=="Eupercaria" | allg_info$last_common_ancestor=="Perciformes"] <- "Percomorphaceae"
## Small sample size in Euacanthomorphacea (62), group with Acanthomorphata (second smallest category) - so that now each category has >100 observations
allg_info$last_common_ancestor[allg_info$last_common_ancestor=="Euacanthomorphacea"] <- "Acanthomorphata"
levels(allg_info$last_common_ancestor) <- c(levels(allg_info$last_common_ancestor),"Singleton")
allg_info$last_common_ancestor[allg_info$last_common_ancestor==""] <- "Singleton"
allg_info$last_common_ancestor <- factor(allg_info$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae","Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta","Singleton"))
allg_info$CNV.x <- as.factor(allg_info$CNV.x)
allg_info$CNV.x <- factor(allg_info$CNV.x, levels=c(1,0))
##############
dupg_info$Cat <- factor(dupg_info$Cat, levels=c("Young","Medium","Old","Singleton"))
dupg_info$ohno <- factor(dupg_info$ohno, levels=c("SSD","ohnolog","single"))
## LCA -Few Eupercaria and Perciformes (<15), group with Percomorphaceae
dupg_info$last_common_ancestor[dupg_info$last_common_ancestor=="Eupercaria" | dupg_info$last_common_ancestor=="Perciformes"] <- "Percomorphaceae"
dupg_info$last_common_ancestor[dupg_info$last_common_ancestor=="Euacanthomorphacea"] <- "Acanthomorphata"
dupg_info$last_common_ancestor <- factor(dupg_info$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae", "Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta",""))
dupg_info$CNVs <- as.factor(dupg_info$CNVs)
dupg_info$CNVs <- factor(dupg_info$CNVs, levels=c(2,1,0))

##########################################################################################
##Filter: only keep considering pairs with high coverage
threshold=10 ##at least 10 sites counted in promoters for each duplicate

a<-allgenes[[1]]
d<-data[[1]]

filt<-data[[1]]
filt$avg_prom_meth<-ifelse((filt$prom_counts>threshold&filt$p_prom_counts>threshold),filt$avg_prom_meth,NA)
filt$div_prom_meth<-ifelse((filt$prom_counts>threshold&filt$p_prom_counts>threshold),filt$div_prom_meth,NA)
ap<-as.data.frame(filt$avg_prom_meth)
dp<-as.data.frame(filt$div_prom_meth)
for (f in 2:length(files)){
  ffilt<-data[[f]]
  ffilt$avg_prom_meth<-ifelse((ffilt$prom_counts>threshold&ffilt$p_prom_counts>threshold),ffilt$avg_prom_meth,NA)
  ffilt$div_prom_meth<-ifelse((ffilt$prom_counts>threshold&ffilt$p_prom_counts>threshold),ffilt$div_prom_meth,NA)
  ap<-cbind(ap,ffilt$avg_prom_meth)
  dp<-cbind(dp,ffilt$div_prom_meth)
  
  a<-rbind(a,allgenes[[f]]) ##all data
  d<-rbind(d,data[[f]]) ##all duplicates
  filt<-rbind(filt,ffilt) ##only filtered
}
colnames(ap)<-meta$Sample_Name
colnames(dp)<-meta$Sample_Name

write.csv(a, "Concatenated.gene.methylation.csv",row.names=F)
write.csv(d, "Concatenated.dup.methylation.csv",row.names=F)
write.csv(ap, "Concatenated.gene.methylation.indivs.csv",row.names=F)
write.csv(dp, "Concatenated.dup.methylation.indivs.csv",row.names=F)

##########################################################################################
# heatmaps
salinity <- meta$treatment
salinity<-gsub("21ppt salinity","gray50",salinity)
salinity<-gsub("2ppt salinity","gray88",salinity)

pdf("Methylation_Heatmap_allgenes_promoter.pdf")
heatmap.2(as.matrix(na.omit(plotallg)),col=rev(morecols(50)),dendrogram="column",trace="none", main="Promoter methylation",ColSideColors=salinity,scale="row")
dev.off()
pdf("Methylation_Heatmap_allgenes_genebody.pdf")
heatmap.2(as.matrix(na.omit(plotallgb)),col=rev(morecols(50)),dendrogram="column",trace="none", main="Genebody methylation",ColSideColors=salinity,scale="row")
dev.off()

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
### GLM Promoter methylation ... look at influence of sex and salinity

full_a<-merge(a,allg_info)
full_d<-merge(d,dupg_info)

fit<-lm(formula = promoter_methylation ~ Cat + sex + treatment, data = full_a)
summary(fit)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
##Perform analysis of average across individuals, set up tables

Ameth<-aggregate(.~ensembl_gene_id, data=a[,names(a) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation")], mean)
Ameth_sex<-aggregate(.~ensembl_gene_id+sex, data=a[,names(a) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","sex")], mean)
Ameth_sal<-aggregate(.~ensembl_gene_id+treatment, data=a[,names(a) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","treatment")], mean)
Dmeth<-aggregate(.~ensembl_gene_id, data=d[,names(d) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation")], mean)
Dmeth_sex<-aggregate(.~ensembl_gene_id+sex, data=d[,names(d) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation","sex")], mean)
Dmeth_sal<-aggregate(.~ensembl_gene_id+treatment, data=d[,names(d) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation","treatment")], mean)

Dmeth$avg_prom_meth<-(Dmeth$promoter_methylation+Dmeth$p_promoter_methylation)/2
Dmeth$div_prom_meth<-(Dmeth$promoter_methylation - Dmeth$p_promoter_methylation) / (Dmeth$promoter_methylation + Dmeth$p_promoter_methylation)
Dmeth$avg_prom_counts<-(Dmeth$prom_counts+Dmeth$p_prom_counts)/2
Dmeth$avg_gb_meth<-(Dmeth$genebody_methylation+Dmeth$p_genebody_methylation)/2
Dmeth$div_gb_meth<-(Dmeth$genebody_methylation - Dmeth$p_genebody_methylation) / (Dmeth$genebody_methylation + Dmeth$p_genebody_methylation)

Dmeth_sex$avg_prom_meth<-(Dmeth_sex$promoter_methylation+Dmeth_sex$p_promoter_methylation)/2
Dmeth_sex$div_prom_meth<-(Dmeth_sex$promoter_methylation - Dmeth_sex$p_promoter_methylation) / (Dmeth_sex$promoter_methylation + Dmeth_sex$p_promoter_methylation)
Dmeth_sex$avg_prom_counts<-(Dmeth_sex$prom_counts+Dmeth_sex$p_prom_counts)/2
Dmeth_sex$avg_gb_meth<-(Dmeth_sex$genebody_methylation+Dmeth_sex$p_genebody_methylation)/2
Dmeth_sex$div_gb_meth<-(Dmeth_sex$genebody_methylation - Dmeth_sex$p_genebody_methylation) / (Dmeth_sex$genebody_methylation + Dmeth_sex$p_genebody_methylation)

Dmeth_sal$avg_prom_meth<-(Dmeth_sal$promoter_methylation+Dmeth_sal$p_promoter_methylation)/2
Dmeth_sal$div_prom_meth<-(Dmeth_sal$promoter_methylation - Dmeth_sal$p_promoter_methylation) / (Dmeth_sal$promoter_methylation + Dmeth_sal$p_promoter_methylation)
Dmeth_sal$avg_prom_counts<-(Dmeth_sal$prom_counts+Dmeth_sal$p_prom_counts)/2
Dmeth_sal$avg_gb_meth<-(Dmeth_sal$genebody_methylation+Dmeth_sal$p_genebody_methylation)/2
Dmeth_sal$div_gb_meth<-(Dmeth_sal$genebody_methylation - Dmeth_sal$p_genebody_methylation) / (Dmeth_sal$genebody_methylation + Dmeth_sal$p_genebody_methylation)

##filter
threshold=10 ##at least 10 sites counted in promoters for each gene
afilter<-a
afilter$promoter_methylation<-ifelse(afilter$prom_counts>=threshold,afilter$promoter_methylation,NA)
dfilter<-d
dfilter$promoter_methylation<-ifelse(dfilter$prom_counts>=threshold,dfilter$promoter_methylation,NA)
dfilter$p_promoter_methylation<-ifelse(dfilter$p_prom_counts>=threshold,dfilter$p_promoter_methylation,NA)

Amethf<-aggregate(.~ensembl_gene_id, data=afilter[,names(afilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation")], mean)
Amethf_sex<-aggregate(.~ensembl_gene_id+sex, data=afilter[,names(afilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","sex")], mean)
Amethf_sal<-aggregate(.~ensembl_gene_id+treatment, data=afilter[,names(afilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","treatment")], mean)
Dmethf<-aggregate(.~ensembl_gene_id, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","p_promoter_methylation","p_prom_counts")], mean)
Dmethf_sex<-aggregate(.~ensembl_gene_id+sex, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","p_promoter_methylation","p_prom_counts","sex")], mean)
Dmethf_sal<-aggregate(.~ensembl_gene_id+treatment, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","p_promoter_methylation","p_prom_counts","treatment")], mean)
Dmethf<-aggregate(.~ensembl_gene_id, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation")], mean)
Dmethf_sex<-aggregate(.~ensembl_gene_id+sex, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation","sex")], mean)
Dmethf_sal<-aggregate(.~ensembl_gene_id+treatment, data=dfilter[,names(dfilter) %in% c("ensembl_gene_id","promoter_methylation","prom_counts","genebody_methylation","p_promoter_methylation","p_prom_counts","p_genebody_methylation","treatment")], mean)

Dmethf$avg_prom_meth<-(Dmethf$promoter_methylation+Dmethf$p_promoter_methylation)/2
Dmethf$div_prom_meth<-(Dmethf$promoter_methylation - Dmethf$p_promoter_methylation) / (Dmethf$promoter_methylation + Dmethf$p_promoter_methylation)
Dmethf$avg_prom_counts<-(Dmethf$prom_counts+Dmethf$p_prom_counts)/2
Dmethf$avg_gb_meth<-(Dmethf$genebody_methylation+Dmethf$p_genebody_methylation)/2
Dmethf$div_gb_meth<-(Dmethf$genebody_methylation - Dmethf$p_genebody_methylation) / (Dmethf$genebody_methylation + Dmethf$p_genebody_methylation)

Dmethf_sex$avg_prom_meth<-(Dmethf_sex$promoter_methylation+Dmethf_sex$p_promoter_methylation)/2
Dmethf_sex$div_prom_meth<-(Dmethf_sex$promoter_methylation - Dmethf_sex$p_promoter_methylation) / (Dmethf_sex$promoter_methylation + Dmethf_sex$p_promoter_methylation)
Dmethf_sex$avg_prom_counts<-(Dmethf_sex$prom_counts+Dmethf_sex$p_prom_counts)/2
Dmethf_sex$avg_gb_meth<-(Dmethf_sex$genebody_methylation+Dmethf_sex$p_genebody_methylation)/2
Dmethf_sex$div_gb_meth<-(Dmethf_sex$genebody_methylation - Dmethf_sex$p_genebody_methylation) / (Dmethf_sex$genebody_methylation + Dmethf_sex$p_genebody_methylation)

Dmethf_sal$avg_prom_meth<-(Dmethf_sal$promoter_methylation+Dmethf_sal$p_promoter_methylation)/2
Dmethf_sal$div_prom_meth<-(Dmethf_sal$promoter_methylation - Dmethf_sal$p_promoter_methylation) / (Dmethf_sal$promoter_methylation + Dmethf_sal$p_promoter_methylation)
Dmethf_sal$avg_prom_counts<-(Dmethf_sal$prom_counts+Dmethf_sal$p_prom_counts)/2
Dmethf_sal$avg_gb_meth<-(Dmethf_sal$genebody_methylation+Dmethf_sal$p_genebody_methylation)/2
Dmethf_sal$div_gb_meth<-(Dmethf_sal$genebody_methylation - Dmethf_sal$p_genebody_methylation) / (Dmethf_sal$genebody_methylation + Dmethf_sal$p_genebody_methylation)

###########################################################################################################################################################################################
##Add info to dataframes
Ameth<-merge(Ameth,allg_info)
Ameth_sex<-merge(Ameth_sex,allg_info)
Ameth_sal<-merge(Ameth_sal,allg_info)
Dmeth<-merge(Dmeth,dupg_info)
Dmeth_sex<-merge(Dmeth_sex,dupg_info)
Dmeth_sal<-merge(Dmeth_sal,dupg_info)
Amethf<-merge(Amethf,allg_info)
Amethf_sex<-merge(Amethf_sex,allg_info)
Amethf_sal<-merge(Amethf_sal,allg_info)
Dmethf<-merge(Dmethf,dupg_info)
Dmethf_sex<-merge(Dmethf_sex,dupg_info)
Dmethf_sal<-merge(Dmethf_sal,dupg_info)

write.csv(Dmethf, "Merged.methylation.filtered.dups2.csv",row.names=F)
write.csv(Amethf, "Merged.methylation.filtered.all2.csv",row.names=F)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
#########plotting
dataprint<-"Combined"

RS.fil<-read.csv("Combined.bismark.cov.gz.features.csv",row.names=1)
ggplot(RS.fil, aes(x=region, y=pct)) +
  geom_violin(aes(fill=region), show.legend=FALSE,draw_quantiles=c(0.5)) +
  scale_x_discrete(labels=c("intergenic_nonpromoter"="Intergenic", "introns"="Introns", "exons"="Exons", "promoters"="Promoters")) +
  labs(x="", y="Percent Methylation") +
  theme_light(base_size = 16)

ggplot(Amethf, aes(y=genebody_methylation,x=promoter_methylation)) + 
  geom_point(alpha=0.1) + stat_smooth(method="lm") + 
  ylab("genebody methylation") + xlab("promoter methylation") +
  theme_light(base_size = 16)

###########################################################################################################################################################################################
##Methylation and expression

pm<-ggplot(Amethf, aes(x=promoter_methylation,y=log2(average))) + 
  geom_point(alpha=0.1)  + stat_smooth(method="lm") + 
  xlab("promoter methylation") + ylab("log2 expression") +
  theme_light(base_size = 16)
gb<-ggplot(Amethf, aes(x=genebody_methylation,y=log2(average))) + 
  geom_point(alpha=0.1) + stat_smooth(method="lm") + 
  xlab("genebody methylation") + ylab("log2 expression") +
  theme_light(base_size = 16)

###########################################################################################################################################################################################
## Methylation per genome position
groups<-RS.fil[grep("group",RS.fil$chr),]
scaffs<-RS.fil[grep("scaffold_",RS.fil$chr),]
scaffs$chr<-"scaffold"
##across chromosomes
aggs<-aggregate(pct ~ region + chr, rbind(groups,scaffs), mean)
aggs$chr<-factor(aggs$chr,levels=c("groupI","groupII","groupIII","groupIV","groupV","groupVI","groupVII","groupVIII","groupIX","groupX","groupXI","groupXII","groupXIII","groupXIV","groupXV","groupXVI","groupXVII","groupXVIII","groupXIX","groupXX","groupXXI","scaffold"))
aggs$region<-factor(aggs$region,levels=c("intergenic_nonpromoter","introns","exons","promoters","TSSes"))
aggs$sex<-"all"

ggplot(aggs, aes(x=region, y=pct,fill=region)) +
  facet_grid(~chr, switch='y') +
  geom_bar(stat = "identity") 

###########################################################################################################################################################################################
## Loop through each indiv to see variation across chromosomes

listoffiles<-c("SRR6410472_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410473_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410474_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410475_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410476_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410477_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410478_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410479_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410480_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410481_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410482_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv","SRR6410483_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.features.csv")
for (f in 1:length(listoffiles)){
  data<-listoffiles[f]
  RS.fil.i<-read.csv(data,row.names=1)
  groups<-RS.fil.i[grep("group",RS.fil.i$chr),]
  scaffs<-RS.fil.i[grep("scaffold_",RS.fil.i$chr),]
  scaffs$chr<-"scaffold"
  aggs.i<-aggregate(pct ~ region + chr, rbind(groups,scaffs), mean)
  aggs.i$chr<-factor(aggs.i$chr,levels=c("groupI","groupII","groupIII","groupIV","groupV","groupVI","groupVII","groupVIII","groupIX","groupX","groupXI","groupXII","groupXIII","groupXIV","groupXV","groupXVI","groupXVII","groupXVIII","groupXIX","groupXX","groupXXI","scaffold"))
  aggs.i$sex<-meta[f,6]
  aggs<-rbind(aggs,aggs.i) ##all data
}

aggs$chr<-factor(aggs$chr,levels=c("groupI","groupII","groupIII","groupIV","groupV","groupVI","groupVII","groupVIII","groupIX","groupX","groupXI","groupXII","groupXIII","groupXIV","groupXV","groupXVI","groupXVII","groupXVIII","groupXIX","groupXX","groupXXI","scaffold"))
aggs$region<-factor(aggs$region,levels=c("intergenic_nonpromoter","introns","exons","promoters","TSSes"))

ggplot(aggs[aggs$region!="TSSes",], aes(x=chr, y=pct,fill=chr)) +
  facet_grid(~region, switch='y') +
  geom_boxplot(outlier.size =1) +
  ylab("methylation") + xlab("chromosome") +
  theme_light(base_size = 16) +
  theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.y = element_text(angle = 180))

###########################################################################################################################################################################################
## LCA trends across duplicate gene pairs

fexp<-ggplot(Dmethf, aes(x=last_common_ancestor, y=log2(avg_exp),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("average expression") + xlab("last common ancestor") 
fgb<-ggplot(Dmethf, aes(x=last_common_ancestor, y=(genebody_methylation),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("genebody methylation") + xlab("last common ancestor") 
fmeth<-ggplot(Dmethf, aes(x=last_common_ancestor, y=abs(avg_prom_meth),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation") + xlab("last common ancestor")

###########################################################################################################################################################################################
###Methylation / Expression vs Categories

##add singletons info for plotting
singletons<-Amethf[Amethf$last_common_ancestor=="Singleton",names(Amethf) %in% c("ensembl_gene_id","promoter_methylation","genebody_methylation","last_common_ancestor","ohno","Cat","average")]
Dmethf_sing<-Dmethf[,names(Dmethf) %in% c("ensembl_gene_id","avg_prom_meth","avg_gb_meth","last_common_ancestor","ohno","Cat","avg_exp")]
names(Dmethf_sing)<-names(singletons)
Dmethf_sing<-rbind(Dmethf_sing,singletons)
Dmethf_sing$logavg<-log2(Dmethf_sing$average)
Dmethf_sing$ohno <- factor(Dmethf_sing$ohno, levels=c("SSD","ohnolog","single"))

########### Stats
fmeth<-ggboxplot(Dmethf_sing, x="Cat", y="promoter_methylation",fill="Cat",varwidth=T,notch=T,outlier.shape=NA,
                 ylab="promoter methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Young","Old"),c("Young","Singleton")),na.rm = TRUE) 

fgb<-ggboxplot(Dmethf_sing, x="Cat", y="genebody_methylation",fill="Cat",varwidth=T,notch=T,
               ylab="genebody methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Young","Old"),c("Young","Singleton")),na.rm = TRUE) 

fexp<-ggboxplot(Dmethf_sing, x="Cat", y="logavg",fill="Cat",varwidth=T,notch=T,
                ylab="average expression",xlab="",legend="none") +
	stat_compare_means(comparisons = list(c("Young","Medium"),c("Young","Old"),c("Young","Singleton")),na.rm = TRUE) 

ggarrange(fmeth,fgb,fexp,
          labels = c("A", "B","C"),
          ncol = 1, nrow = 3)

######################################
## with ohno
ggplot(Dmethf_sing[], aes(x=Cat, y=logavg,fill=Cat)) +
  facet_grid(~ohno) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE,varwidth = T) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("average expression") + xlab("last common ancestor") 
ggplot(Dmethf_sing, aes(x=Cat, y=(genebody_methylation),fill=Cat)) +
  facet_grid(~ohno) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE,varwidth = T) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("genebody methylation") + xlab("last common ancestor") +
  #  theme(plot.margin=unit(c(1,1,1,2),"cm")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggplot(Dmethf_sing, aes(x=Cat, y=abs(promoter_methylation),fill=Cat)) +
  facet_grid(~ohno) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE,varwidth = T) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation") + xlab("last common ancestor") +
  #theme(plot.margin=unit(c(1,1,1,2),"cm")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

######################################
## only SSDs
fmeth<-ggboxplot(Dmethf_sing[Dmethf_sing$ohno=="SSD",], x="Cat", y="promoter_methylation",fill="Cat",varwidth=T,notch=T,outlier.shape=NA,
                 ylab="promoter methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

fgb<-ggboxplot(Dmethf_sing[Dmethf_sing$ohno=="SSD",], x="Cat", y="genebody_methylation",fill="Cat",varwidth=T,notch=T,outlier.shape=NA,
               ylab="genebody methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

fexp<-ggboxplot(Dmethf_sing[Dmethf_sing$ohno=="SSD",], x="Cat", y="logavg",fill="Cat",varwidth=T,notch=T,outlier.shape=NA,
                ylab="average expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

ggarrange(fmeth,fgb,fexp,
          labels = c("A", "B","C"),
          ncol = 1, nrow = 3)
ggsave(paste(dataprint,"FigureS5_stats_SSD.pdf",sep="."),height=8,width=4)


###########################################################################################################################################################################################
###per individual
afilter_info<-merge(afilter[,names(afilter) %in% c("ensembl_gene_id","promoter_methylation","genebody_methylation","sample","sex","treatment")],allg_info[,names(allg_info) %in% c("ensembl_gene_id","last_common_ancestor","ohno","Cat","average")])
singletons<-afilter_info[afilter_info$last_common_ancestor=="Singleton",]
dfilter_info<-merge(dfilter[,names(dfilter) %in% c("ensembl_gene_id","avg_prom_meth","avg_gb_meth","sample","sex","treatment")],dupg_info[,names(dupg_info) %in% c("ensembl_gene_id","last_common_ancestor","ohno","Cat","avg_exp")])
names(dfilter_info)<-names(afilter_info)
dfilter_info<-rbind(dfilter_info,singletons)

fpm<-ggplot(dfilter_info, aes(x=sample, y=promoter_methylation,fill=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T,show.legend=TRUE,outlier.shape=NA) +
  theme_light(base_size = 14) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("") + ylab("promoter methylation")

fgb<-ggplot(dfilter_info,  aes(x=sample, y=genebody_methylation,fill=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T,show.legend=TRUE,outlier.shape=NA) +
  theme_light(base_size = 14) +
  xlab("") + ylab("genebody methylation") 

ggarrange(fpm,fgb,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

###########################################################################################################################################################################################
## dS, remove very large dS (>3)
dS<-Dmethf[Dmethf$dS<3.01 & !is.na(Dmethf$dS),]
dS$bin2<-bin(dS$dS,nbins=10,method="content") ##create equal-sized bins (method=content)

##dS Separate into bins!! redo bins where all genes also have methylation data
dSs<-dS[!is.na(dS$avg_prom_meth),]
dSs$bin2<-bin(dSs$dS,nbins=10,method="content") ##create bins,equal-sized with method="content"

ggplot(dS, aes(x=last_common_ancestor, y=dS,fill=last_common_ancestor)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 14)  +  ylab("dS") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12)) +
  xlab("dS") +
  theme(legend.position="none")

fexp<-ggplot(dS, aes(x=bin2, y=log2(avg_exp),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 14)  +  ylab("average expression") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12)) +
  xlab("dS") +
  theme(legend.position="none")
fgb<-ggplot(dS, aes(x=bin2, y=abs(avg_gb_meth),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 14)  +  ylab("genebody methylation") + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
fmeth<-ggplot(dS, aes(x=bin2, y=abs(avg_prom_meth),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 14) +  ylab("promoter methylation") + 
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggarrange(fmeth,fgb,fexp,
          labels = c("A", "B","C"),
          ncol = 1, nrow = 3, heights=c(1,1,1.67))

###########################################################################################################################################################################################
### Compare with sig differentiated salinity from Metzger - differentially methylated cytosines
#DMCs salinity in promoters : 71/(1186+71) = 5.6%, similar to proportion of sites in genome:  24,921,487/446,627,861
#DMCs sex in promoters : 1052/(18563) , similar to proportion of sites in genome:  24,921,487/446,627,861
#DMCs sex in promoters: 602 genes

DMC<-read.table("DMCs_Metzger_Salinity_relabeledchroms.bed.Promoters.bed") ##calculated DMCs in promoters in bed format, using bedtools
names(DMC)<-c("chr","start","end","pvalue","qvalue","meth.diff","pchr","pstart","pend","ensembl_gene_id","tid","strand","gstart","gend","biotype","CG","CGoverlap","length","CGpercentage","overlap")
##promoters with DMCs
DMCprom<-as.data.frame(table(DMC$ensembl_gene_id))
names(DMCprom)[1]<-"ensembl_gene_id"
DMCpromg<-DMCprom[-1,]
DMCtypes<-merge(Amethf,DMCpromg,all.x=T)

##DMCs tend to be in promoters of old duplicate genes
Props<-rbind(genome=table(DMCtypes$Cat)/sum(table(DMCtypes$Cat)),DMCs=table(DMCtypes$Cat[!is.na(DMCtypes$Freq)])/sum(table(DMCtypes$Cat[!is.na(DMCtypes$Freq)])))
dd<-melt(Props)
names(dd)<-c("Region","Cat","Proportion")
ggplot(dd, aes(x=Cat,y=Proportion,fill=Region)) +
  geom_bar(stat="identity",position=position_dodge()) +
  theme_light(base_size = 16)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
##Redo analyses without DMCs, but so few DMCs ... doesn't affect results
ggboxplot(Dmethf_sing[!(Dmethf_sing$ensembl_gene_id %in% DMCpromg$ensembl_gene_id),], x="Cat", y="promoter_methylation",fill="Cat",varwidth=T,notch=T,outlier.shape=NA,
                 ylab="promoter methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Young","Old"),c("Young","Singleton")),na.rm = TRUE)
ggboxplot(Dmethf_sing[!(Dmethf_sing$ensembl_gene_id %in% DMCpromg$ensembl_gene_id),], x="Cat", y="logavg",fill="Cat",varwidth=T,notch=T,
                ylab="average expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Young","Old"),c("Young","Singleton")),na.rm = TRUE) 
ggplot(dS[!(dS$ensembl_gene_id %in% DMCpromg$ensembl_gene_id | dS$paralogue_ensembl_gene_id %in% DMCpromg$ensembl_gene_id),], aes(x=bin2, y=abs(avg_prom_meth),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 14) +  ylab("promoter methylation") + 
  theme(legend.position="none") 
ggplot(dfilter_info[!(dfilter_info$ensembl_gene_id %in% DMCpromg$ensembl_gene_id),], aes(x=sample, y=promoter_methylation,fill=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T,show.legend=TRUE,outlier.shape=NA) +
  ylab("promoter methylation")
###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

#Divergence between duplicates

fexp<-ggplot(Dmethf, aes(x=last_common_ancestor, y=abs(div_exp),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25, notch=T,show.legend=FALSE,outlier.shape=NA) +
  theme_light(base_size = 16) + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +
  xlab("last common ancestor") + ylab("expression divergence") 
fmeth<-ggplot(Dmethf, aes(x=last_common_ancestor, y=abs(div_prom_meth),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25, notch=T,show.legend=FALSE,outlier.shape=NA) +
  theme_light(base_size = 16) + 
  xlab("last common ancestor") + ylab("methylation divergence") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#############################
fexpdS<-ggplot(dS, aes(x=bin2, y=abs(div_exp),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 16) +
  xlab("dS") + ylab("expression divergence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +
  theme(legend.position="none")
fmethdS<-ggplot(dS, aes(x=bin2, y=abs(div_prom_meth),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  theme_light(base_size = 16) +
  xlab("dS") + ylab("methylation divergence") +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggarrange(fmeth,fmethdS,fexp,fexpdS,
          labels = c("A", "C","B","D"),
          ncol = 2, nrow = 2, heights=c(1,1.67))
ggsave(paste(dataprint,"Figure3.pdf",sep="."),width=8)

###########################################################################################################################################################################################
# Divergence among Categories
Dmethf_div<-Dmethf
Dmethf_div$abs_div_meth<-abs(Dmethf_div$div_prom_meth)
Dmethf_div$abs_div_exp<-abs(Dmethf_div$div_exp)

fmeth<-ggboxplot(Dmethf_div, x="Cat", y="abs_div_meth",fill="Cat",varwidth=F,notch=T,
          ylab="divergence methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
fexp<-ggboxplot(Dmethf_div, x="Cat", y="abs_div_exp",fill="Cat",varwidth=F,notch=T,
          ylab="divergence expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

ggarrange(fmeth,fexp,
          labels = c("A", "B"),
          ncol = 1, nrow = 2, heights=c(1,1,1.67))

###only SSDs
fmeth<-ggboxplot(Dmethf_div[Dmethf_div$ohno=="SSD",], x="Cat", y="abs_div_meth",fill="Cat",varwidth=F,notch=T,
                 ylab="divergence methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
fexp<-ggboxplot(Dmethf_div[Dmethf_div$ohno=="SSD",], x="Cat", y="abs_div_exp",fill="Cat",varwidth=F,notch=T,
                ylab="divergence expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

ggarrange(fmeth,fexp,
          labels = c("A", "B"),
          ncol = 1, nrow = 2, heights=c(1,1,1.67))

###########################################################################################################################################################################################
###Every individual - divergence in promoter methylation
#add gene info
Dmethf_div_ind<-cbind(dupg_info[,c(1,3,8,9,12)],dp)
Dmethf_div_ind<-melt(Dmethf_div_ind)
Dmethf_div_ind$variable<-factor(Dmethf_div_ind$variable,levels=sort(levels(Dmethf_div_ind$variable)))

ggplot(Dmethf_div_ind, aes(x=variable, y=abs(value),fill=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T,show.legend=TRUE,outlier.shape=NA) +
  theme_light(base_size = 14) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation divergence") 

###remove low dS
Dmethf_div_ind<-cbind(dupg_info[,c(1,3,5,8,9,12)],dp)
Dmethf_div_ind<-Dmethf_div_ind[Dmethf_div_ind$dS>0.05,]
Dmethf_div_ind<-melt(Dmethf_div_ind)
Dmethf_div_ind<-Dmethf_div_ind[Dmethf_div_ind$variable!="dS",]
Dmethf_div_ind$variable<-factor(Dmethf_div_ind$variable,levels=sort(levels(Dmethf_div_ind$variable)))

ggplot(Dmethf_div_ind, aes(x=variable, y=abs(value),fill=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T,show.legend=TRUE,outlier.shape=NA) +
  theme_light(base_size = 14) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation divergence") 

###########################################################################################################################################################################################
#correlation of divergence in meth and exp
ggplot(Dmethf, aes(x=div_prom_meth,y=div_exp)) + 
  #geom_point(alpha=0.1) + stat_smooth(method="loess") + 
  geom_point(alpha=0.1) + stat_smooth(method="lm") + 
  xlab("divergence in promoter methylation") + ylab("divergence in expression") +
  theme_light(base_size = 16)
cor.test(Dmethf$div_prom_meth, Dmethf$div_exp,na.rm=TRUE)

#########Correct for expression
### includ others in a glm?
fit<-lm(formula = div_prom_meth ~ Cat + ohnolog + dS + promoter_methylation, data = Dmethf)
summary(fit)

###########################################################################################################################################################################################
# Tissue expression
avgs_all$avg<-(avgs_all$average+avgs_all$p_average)/2
avgs_all$avg_em<-(avgs_all$em+avgs_all$p_em)/2
avgs_all$avg_lm<-(avgs_all$lm+avgs_all$p_lm)/2
avgs_all$avg_mg<-(avgs_all$mg+avgs_all$p_mg)/2
avgs_all$avg_eo<-(avgs_all$eo+avgs_all$p_eo)/2
avgs_all$avg_24<-(avgs_all$X24+avgs_all$p_X24)/2
avgs_all$avg_HK<-(avgs_all$HK+avgs_all$p_HK)/2
avgs_all$avg_SP<-(avgs_all$SP+avgs_all$p_SP)/2
avgs_all$avg_ovaries<-(avgs_all$ovaries+avgs_all$p_ovaries)/2
avgs_all$avg_testes<-(avgs_all$testes+avgs_all$p_testes)/2

avgs_all$div<-(avgs_all$average-avgs_all$p_average)/(avgs_all$average+avgs_all$p_average)
avgs_all$div_em<-(avgs_all$em-avgs_all$p_em)/(avgs_all$em+avgs_all$p_em)
avgs_all$div_lm<-(avgs_all$lm-avgs_all$p_lm)/(avgs_all$lm+avgs_all$p_lm)
avgs_all$div_mg<-(avgs_all$mg-avgs_all$p_mg)/(avgs_all$mg+avgs_all$p_mg)
avgs_all$div_eo<-(avgs_all$eo-avgs_all$p_eo)/(avgs_all$eo+avgs_all$p_eo)
avgs_all$div_24<-(avgs_all$X24-avgs_all$p_X24)/(avgs_all$X24+avgs_all$p_X24)
avgs_all$div_HK<-(avgs_all$HK-avgs_all$p_HK)/(avgs_all$HK+avgs_all$p_HK)
avgs_all$div_SP<-(avgs_all$SP-avgs_all$p_SP)/(avgs_all$SP+avgs_all$p_SP)
avgs_all$div_ovaries<-(avgs_all$ovaries-avgs_all$p_ovaries)/(avgs_all$ovaries+avgs_all$p_ovaries)
avgs_all$div_testes<-(avgs_all$testes-avgs_all$p_testes)/(avgs_all$testes+avgs_all$p_testes)

avgs<-melt(avgs_all)
avgs$Cat <- factor(avgs$Cat, levels=c("Young","Medium","Old"))
avgs$last_common_ancestor <- factor(avgs$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae","Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta","Singleton"))

##expression across tissues by age
e1<-ggplot(avgs[grep("^avg",avgs$variable),], aes(x=Cat, y=log2(value),col=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T) +
  labs(title="", x="", y="log2 expression") +
  theme_light(base_size = 16) + theme(legend.position="top") +
  ylab("expression") +
  facet_grid(~variable) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
e2<-ggplot(avgs[grep("^div",avgs$variable),], aes(x=Cat, y=abs(value),col=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T) +
  labs(title="", x="", y="log2 expression") +
  theme_light(base_size = 16) + theme(legend.position="none") +
  ylab("expression divergence") +
  facet_grid(~variable) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggarrange(e1,e2,
          labels = c("A","B"),
          ncol = 1, nrow = 2,heights=c(1.35,1.0))

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################
##CpG and GC expectations
##Get gene and transcript ids
gg<-read.table("Gasterosteus_aculeatus.Genes.bed",header=T,comment.char="",sep="\t")
colnames(gg)[c(1:5,9)] <- c("seqname", "start", "end", "ensembl_gene_id","ensembl_transcript_id","biotype")

##Get GC from bedtools: 
#bedtools nuc -fi Gasterosteus_aculeatus.dna_sm.toplevel.fa -bed Gasterosteus_aculeatus.Genes.first.bed.promoters.adjusted > Gasterosteus_aculeatus.Genes.first.bed.promoters.adjusted.GC
GCprom<-read.table("Gasterosteus_aculeatus.gff3.genePred.bed.filtered.promoters.fixed.GC",header=T,comment.char="")
names(GCprom)<-c("chrom","start","end","len","strand","x","ensembl_transcript_id","pct_AT","pct_GC","A","C","G","T")

#sed 's/ /	/g'  CG.stickleback.bed | intersectBed -b stdin -a Gasterosteus_aculeatus.gff3.genePred.bed.filtered.promoters.fixed -c > Gasterosteus_aculeatus.gff3.genePred.bed.filtered.promoters.fixed.CpGcount
CGprom<-read.table("Gasterosteus_aculeatus.gff3.genePred.bed.filtered.promoters.fixed.CpGcount",header=F,comment.char="")
names(CGprom)<-c("chrom","start","end","len","strand","x","ensembl_transcript_id","CpG")
CGprom$pct_CpG<-CGprom$CpG/2001
############

CpGmeth<-merge(CGprom,GCprom[,c(1:3,7:13)])
CpGmeth<-merge(CpGmeth,gg[,c(4:5)])
CpGmeth$CGexp<-(CpGmeth$C*CpGmeth$G)/(CpGmeth$len-1)
CpGmeth$CGratio<-CpGmeth$CpG/CpGmeth$CGexp
CGmeth<-merge(CpGmeth,Amethf)

pm<-ggplot(CGmeth[!is.na(CGmeth$promoter_methylation),], aes(x=promoter_methylation, y=log2(CGratio))) +
  geom_point(alpha=0.4,size=0.5) +
  stat_smooth(method="lm") +
  ylab("log CpG ratio") + xlab("promoter methylation") +
  theme_light(base_size = 16)
gb<-ggplot(CGmeth[!is.na(CGmeth$genebody_methylation),], aes(x=genebody_methylation, y=log2(CGratio))) +
  geom_point(alpha=0.4,size=0.5) +
  stat_smooth(method="lm") +
  ylab("log CpG ratio") + xlab("genebody methylation") +
  theme_light(base_size = 16)
ex<-ggplot(CGmeth[!is.na(CGmeth$average),], aes(x=log2(average), y=log2(CGratio))) +
  geom_point(alpha=0.4,size=0.5) +
  stat_smooth(method="lm") +
  ylab("log CpG ratio") + xlab("log expression") +
  theme_light(base_size = 16)

############################################################################################################
####Duplicate genes
###Add to duplicates list, first including singletons
CGmeth_f<-merge(Dmethf_sing,CpGmeth)
##add CNVs
CGmeth_f<-merge(CGmeth_f,Dmethf[,names(Dmethf) %in% c("ensembl_gene_id","CNVs")],all.x=T)
CGmeth_f$CNVs[is.na(CGmeth_f$CNVs)]<-0

###CpG% in promoters (CpG/length)
cg1<-ggboxplot(CGmeth_f, x="Cat", y="pct_CpG",fill="Cat",varwidth=F,notch=T,ylab="CpG percentage",xlab="",legend="none",outlier.size=0.5) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
cg2<-ggboxplot(CGmeth_f, x="Cat", y="pct_GC",fill="Cat",varwidth=F,notch=T,ylab="percent GC",xlab="",legend="none",outlier.size=0.5) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
CGmeth_flog<-CGmeth_f
CGmeth_flog$CGratio<-log2(CGmeth_flog$CGratio)
cg3<-ggboxplot(CGmeth_flog, x="Cat", y="CGratio",fill="Cat",varwidth=F,notch=T,ylab="log CpG Obs/Exp ratio",xlab="",legend="none",outlier.size=0.5) +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 

ggarrange(cg1,cg2,cg3,
          labels = c("A", "B","C"),
          ncol = 3, nrow = 1)

fit<-lm(promoter_methylation ~ CGratio + logavg, data=CGmeth_f[is.finite(CGmeth_f$logavg),])
summary(fit)

##########################################################################################
##add paralog pair, no singletons
CGmeth_d<-merge(CpGmeth,Dmethf)
CGmeth_d<-merge(CGmeth_d,allg_info)
dupCNV<-allg_info[,c("ensembl_gene_id","paralogue_ensembl_gene_id","CNV.x")]
names(dupCNV)<-c("paralogue_ensembl_gene_id","ensembl_gene_id","CNV.y")
dupCNV<-unique(dupCNV)
CGmeth_d<-merge(CGmeth_d,dupCNV) ##add paralog CNV
dup<-CpGmeth[,names(CpGmeth) %in% c("CpG","pct_CpG","pct_GC","ensembl_gene_id","CGexp","CGratio")]
names(dup)<-paste("p_",names(dup),sep="")
names(dup)[4]<-"paralogue_ensembl_gene_id"
CGmeth_d<-merge(CGmeth_d,dup) ##add paralog CG/GC stuff

###Add average and divergence
CGmeth_d$avg_pct_CpG<-(CGmeth_d$pct_CpG+CGmeth_d$p_pct_CpG)/2
CGmeth_d$avg_pct_GC<-(CGmeth_d$pct_GC+CGmeth_d$p_pct_GC)/2
CGmeth_d$avg_pct_CGratio<-(CGmeth_d$CGratio+CGmeth_d$p_CGratio)/2

CGmeth_d$div_pct_CpG<-(CGmeth_d$pct_CpG-CGmeth_d$p_pct_CpG)/(CGmeth_d$pct_CpG+CGmeth_d$p_pct_CpG) ###add divergence
CGmeth_d$div_pct_GC<-(CGmeth_d$pct_GC-CGmeth_d$p_pct_GC)/(CGmeth_d$pct_GC+CGmeth_d$p_pct_GC) ###add divergence
CGmeth_d$div_pct_CGratio<-(CGmeth_d$CGratio-CGmeth_d$p_CGratio)/(CGmeth_d$CGratio+CGmeth_d$p_CGratio) ###add divergence

CGmeth_d$abs_div_pct_CpG<-abs(CGmeth_d$div_pct_CpG)
CGmeth_d$abs_div_pct_GC<-abs(CGmeth_d$div_pct_GC)
CGmeth_d$abs_div_pct_CGratio<-abs(CGmeth_d$div_pct_CGratio)

### Avg by cat
cg1<-ggboxplot(CGmeth_d, x="Cat", y="avg_pct_CpG",fill="Cat",varwidth=F,notch=T,ylab="avg CpG promoter",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
cg2<-ggboxplot(CGmeth_d, x="Cat", y="avg_pct_GC",fill="Cat",varwidth=F,notch=T,ylab="avg pct GC",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
cg3<-ggboxplot(CGmeth_d, x="Cat", y="avg_pct_CGratio",fill="Cat",varwidth=F,notch=T,ylab="avg CG ratio",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
### Div by cat
cg4<-ggboxplot(CGmeth_d, x="Cat", y="abs_div_pct_CpG",fill="Cat",varwidth=F,notch=T,ylab="div CpG promoter",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
cg5<-ggboxplot(CGmeth_d, x="Cat", y="abs_div_pct_GC",fill="Cat",varwidth=F,notch=T,ylab="div pct GC",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
cg6<-ggboxplot(CGmeth_d, x="Cat", y="abs_div_pct_CGratio",fill="Cat",varwidth=F,notch=T,ylab="div CG ratio",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium"),c("Medium","Old"),c("Young","Old")),na.rm = TRUE) 
ggarrange(cg1,cg2,cg3,cg4,cg5,cg6,
          labels = c("A", "B","C","D","E","F"),
          ncol = 3, nrow = 2)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
# CNVs and methylation / expression

Amethf$log<-log2(Amethf$average+0.001)

ee<-Amethf[Amethf$Cat=="Singleton",]
p1<-as.character(Dmethf$ensembl_gene_id)
p2<-as.character(Dmethf$paralogue_ensembl_gene_id)
ee<-rbind(ee,Amethf[Amethf$ensembl_gene_id %in% c(p1,p2),])

ee<-melt(Amethf[,colnames(Amethf) %in% c("average","CNV.x","ohno","Cat")])
ee$value<-log2(ee$value)
a1<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",varwidth=F,notch=T,outlier.size=.25,ylab="expression",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(-11,19)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE) 
a1_full<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="ohno",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="expression",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(-11,19)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
a1_cat<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="Cat",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="expression",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(-11,19)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
ee$ohno <- factor(ee$ohno, levels=c("SSD","ohnolog","single"))
a1_full_1<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="ohno",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="expression",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(-11,19)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)

ee<-melt(Amethf[,colnames(Amethf) %in% c("promoter_methylation","CNV.x","ohno","Cat")])
ee$value<-abs(ee$value)
a2<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(0,105)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
a2_full<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="ohno",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(0,105)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
a2_cat<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="Cat",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(0,105)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by=c("Cat","ohno"),nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(0,105)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)
ee$ohno <- factor(ee$ohno, levels=c("SSD","ohnolog","single"))
a2_full_1<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",facet.by="ohno",nrow=1,varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue"),ylim=c(0,105)) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE)

ee<-melt(Amethf[,colnames(Amethf) %in% c("genebody_methylation","CNV.x")])
ee$value<-abs(ee$value)
a3<-ggboxplot(ee, x="CNV.x", y="value",fill="CNV.x",varwidth=F,notch=T,outlier.size=.25,ylab="genebody methylation",xlab="",legend="none",palette = c("#E7B800","cornflowerblue")) +
  stat_compare_means(comparisons = list(c("0","1")),na.rm = TRUE) 

ggarrange(a2,g1,g3,a1,g2,g4,
          labels = c("A","C","E","B","D","F"),
          ncol = 3, nrow = 2)

ggarrange(a2_full,a2_cat,a1_full,a1_cat,
          labels = c("A","C","B","D"),
          ncol = 2, nrow = 2)

ggarrange(a2_full_1,a1_full_1,
          labels = c("A","B"),
          ncol = 1, nrow = 2)

fit<-lm(formula = promoter_methylation ~ CNV.x + log, data = Amethf)
summary(fit)

###########################################################################################################################################################################################
ggplot(Dmethf, aes(x=CNVs, y=promoter_methylation,fill=Cat)) +
  facet_grid(ohno~Cat) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE,varwidth = F) +
  theme_light(base_size = 16) 

ggplot(Amethf, aes(x=CNV.x, y=promoter_methylation,fill=Cat)) +
  facet_grid(ohno~Cat) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE,varwidth = F) +
  theme_light(base_size = 16) 

##divergence
p1<-as.character(Dmethf$ensembl_gene_id)
p2<-as.character(Dmethf$paralogue_ensembl_gene_id)
ddd<-Amethf[Amethf$ensembl_gene_id %in% c(p1,p2),]
ggboxplot(ddd, x="CNV.x", y="log",facet.by=c("ohnolog","Cat"),fill="Cat",varwidth=F,notch=T,outlier.size=.25,ylab="log2 expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("1","0")),na.rm = TRUE,vjust=0.5) 
ggboxplot(ddd, x="CNV.x", y="promoter_methylation",facet.by=c("ohnolog","Cat"),fill="Cat",varwidth=F,notch=T,outlier.size=.25,ylab="promoter methylation",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("1","0")),na.rm = TRUE,vjust=0.5) 

