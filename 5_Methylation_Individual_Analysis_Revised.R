##########################################################################################
#R script to determine promoter-wide methylation of each gene and gene-wide methylation for an individual
#Needs genome and transcription
#usage: Rscript Methylation_Individual_Analysis_Revised.R InputFile.cov.gz

##Assumes already created duplicate gene list and TranscriptomeBigScreen.R to generate expression data from big screen data
#Incorporates the MethylationDuplicateAnalysis.R at the end

library(GenomicRanges)
library(genomation)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggpubr) ##for pub-quality figs
library(OneR) ##for binning
source("ggplot_smooth_func.R") ##adds regression equation to plot
library(data.table) ##allows tstrsplit

##########################################################################################
# Annotations
# all genes with ortho/para, CNV, etc. and expression data, as well as paralog pairs
genes<-read.csv("Transcriptome.TPM.dup.div.exp.csv",row.names=1) ###Info on stickleback gene expression and duplicate gene expression divergence
ppairs<-read.csv("SticklebackGeneParalogPairs.v99.csv",row.names=1) ###Orth/Para info for each gene pair combined without repeats
bed<-read.csv("Gasterosteus_aculeatus.Genes.bed.csv",row.names=1) ###Gene information
### Only keep protein coding
genes<-genes[genes$Gene=="protein_coding",]
ppairs<-merge(ppairs,genes)
bed<-bed[bed$biotype=="protein_coding",]
##########################################################################################
# Granges features
bed.gr<-as(bed,"GRanges")

# Features: promoter = 1500bp upstream, 500bp downstream
features.gr <- readTranscriptFeatures("Gasterosteus_aculeatus.gff3.genePred.bed.filtered",
                                   up.flank=1501, down.flank=500, unique.prom=FALSE,remove.unusual=FALSE) ##remove.unusual F to keep scaffolds in ENSEMBL version
write.table(features.gr$promoters,"Gasterosteus_aculeatus.gff3.genePred.bed.filtered.promoters",row.names=F,quote=F,sep="\t",col.names=F)

###########################################################################################
###########################################################################################
###########################################################################################

##Load the data from bismark
data<-"Sample.bismark.cov.gz" 
RS <- read.table(gzfile(data), col.names = c("chr", "start", "end", "pct", "numCs", "numTs"))
# create coverage column and filter by coverage
RS$coverage <- RS$numCs + RS$numTs
RS.fil <- RS[RS$coverage >= 10, ] ##filter by coverage
RS.fil <- RS.fil[RS.fil$coverage <= quantile(RS$coverage, .999), ] # 99.9th percentile calculation
RS.gr <- as(RS.fil, "GRanges") # GRanges object needed for findOverlap()

##Add annotations to methylation dataset
geneOverlaps <- as.data.frame(findOverlaps(RS.gr, bed.gr))

RS.fil[geneOverlaps$queryHits, "ensembl_gene_id"] <- as.character(bed.gr[geneOverlaps$subjectHits,]$ensembl_gene_id)
RS.fil[geneOverlaps$queryHits, "ensembl_transcript_id"] <- as.character(bed.gr[geneOverlaps$subjectHits,]$ensembl_transcript_id)
RS.fil$region <- factor("intergenic_nonpromoter", levels=c("intergenic_nonpromoter", "TSSes", "introns", "exons", "promoters")) # default region; later, sites in other regions will have this column overwritten
RS.fil[ , c("feature_start", "feature_end", "feature_score")] <- NA

##loop to annotate, precedence of promoters > exons > TSSes > introns > intergenic
feats <- c("introns", "TSSes", "exons", "promoters")
for (i in 1:length(feats)) {
  feature.gr <- as(features.gr[[feats[i]]], "GRanges")
  feature.df <- as.data.frame(feature.gr)# ensembl_transcript_id, ensembl_gene_id
  colnames(feature.df)[7] <- "ensembl_transcript_id"
  feature.df <- join(feature.df, bed[,4:5])
  
  overlaps <- as.data.frame(findOverlaps(RS.gr, feature.gr))
  RS.fil[overlaps$queryHits, "region"] <- feats[i]
  RS.fil[overlaps$queryHits, "feature_start"] <- start(ranges(feature.gr[overlaps$subjectHits, ]))
  RS.fil[overlaps$queryHits, "feature_end"] <- end(ranges(feature.gr[overlaps$subjectHits, ]))
  RS.fil[overlaps$queryHits, "ensembl_gene_id"] <- as.character(feature.df[overlaps$subjectHits, "ensembl_gene_id"])
  RS.fil[overlaps$queryHits, "ensembl_transcript_id"] <- feature.df[overlaps$subjectHits, "ensembl_transcript_id"]
  RS.fil[overlaps$queryHits, "feature_score"] <- feature.gr[overlaps$subjectHits, ]$score
}
write.csv(RS.fil, paste(data,"features.csv",sep="."))

##########################################################################################
# make table - promoter region
prom <- aggregate(pct ~ feature_start + ensembl_gene_id, RS.fil[RS.fil$region=="promoters", ], mean)
colnames(prom)[3] <- "promoter_methylation"
counts <- aggregate(pct ~ feature_start + ensembl_gene_id, RS.fil[RS.fil$region=="promoters", ], length) ##number of observations 
prom$prom_counts<-counts$pct
prom[,1] <- NULL

# make table - gene body region
genebody <- aggregate(pct ~ ensembl_gene_id, RS.fil[RS.fil$region=="TSSes" | RS.fil$region=="exons" | RS.fil$region=="introns", ], mean)
colnames(genebody)[2] <- "genebody_methylation"
counts <- aggregate(pct ~ ensembl_gene_id, RS.fil[RS.fil$region=="TSSes" | RS.fil$region=="exons" | RS.fil$region=="introns", ], length)
genebody$genebody_counts<-counts$pct

## Make methylation file for analysis and plotting
meth <- merge(prom, genebody, by="ensembl_gene_id",all=T)
write.csv(meth, paste(data,"meth.prom.gene.csv",sep="."))

###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################

## ALL GENES and ADDING DUPLICATE GENE DATA
meth_p<-meth
names(meth_p)<-paste("p_",names(meth_p),sep="")
names(meth_p)[1]<-"paralogue_ensembl_gene_id"

full<-merge(genes,meth,all.x=T)
full<-merge(full,meth_p,all.x=T)
full<-full[,c(2,1,3:ncol(full))] #reorder first 2 columns

## Divergence in methylation
full$avg_prom_meth <- rowMeans(as.data.frame(full)[c("promoter_methylation", "p_promoter_methylation")])
full$diff_prom_meth <- (full$promoter_methylation - full$p_promoter_methylation)
full$div_prom_meth <- (full$promoter_methylation - full$p_promoter_methylation) / (full$promoter_methylation + full$p_promoter_methylation)
full$sum_prom_counts <- rowSums(as.data.frame(full)[c("prom_counts", "p_prom_counts")])

full$avg_gb_meth <- rowMeans(as.data.frame(full)[c("genebody_methylation", "p_genebody_methylation")])
full$diff_gb_meth <- (full$genebody_methylation - full$p_genebody_methylation)
full$div_gb_meth <- (full$genebody_methylation - full$p_genebody_methylation) / (full$genebody_methylation + full$p_genebody_methylation)

## write out results
full_all<-full
write.csv(full_all, paste(data,"meth.prom.gene.full.csv",sep="."))
full<-merge(ppairs,full_all) ##only keep paralog pairs, each pair occurs only once
write.csv(full, paste(data,"meth.prom.gene.dups.csv",sep="."))

####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
# Order factors for plotting
## CNVs
full_all$CNVs<-as.factor(full_all$CNVs)
full_all$CNVs <- factor(full_all$CNVs, levels=c(2,1,0))
full_all$CNV.x<-as.factor(full_all$CNV.x)
full_all$CNV.x <- factor(full_all$CNV.x, levels=c(1,0))
## Category
full_all$Cat <- factor(full_all$Cat, levels=c("Young","Medium","Old","Singleton"))
full_all$ohno <- factor(full_all$ohno, levels=c("SSD","ohnolog","single"))
## LCA -Few Eupercaria and Perciformes (<15), group with Percomorphaceae
full_all$last_common_ancestor[full_all$last_common_ancestor=="Eupercaria" | full_all$last_common_ancestor=="Perciformes"] <- "Percomorphaceae"
levels(full_all$last_common_ancestor) <- c(levels(full_all$last_common_ancestor),"Singleton")
full_all$last_common_ancestor[full_all$last_common_ancestor==""] <- "Singleton"
full_all$last_common_ancestor <- factor(full_all$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae","Euacanthomorphacea", "Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta","Singleton"))

##############
# Order factors for plotting
## CNVs
full$CNVs<-as.factor(full$CNVs)
full$CNVs <- factor(full$CNVs, levels=c(2,1,0))
full$CNV.x<-as.factor(full$CNV.x)
full$CNV.x <- factor(full$CNV.x, levels=c(1,0))
## Category
full$Cat <- factor(full$Cat, levels=c("Young","Medium","Old"))
full$ohno <- factor(full$ohno, levels=c("SSD","ohnolog"))
## LCA -Few Eupercaria and Perciformes (<15), group with Percomorphaceae
full$last_common_ancestor[full$last_common_ancestor=="Eupercaria" | full$last_common_ancestor=="Perciformes"] <- "Percomorphaceae"
full$last_common_ancestor <- factor(full$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae","Euacanthomorphacea", "Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta",""))

##only consider if promoter count is >= min
mincountthreshold<-10
full<-full[full$prom_counts>=mincountthreshold & full$p_prom_counts>=mincountthreshold,]
full_all<-full_all[full_all$prom_counts>=mincountthreshold & full_all$p_prom_counts>=mincountthreshold,]

##########################################################################################
## dS
##dS Separate into bins!! redo bins where all genes also have methylation data

##ignore dS>3
dS<-full[full$dS<3.01 & !is.na(full$dS),]
dS$bin2<-bin(dS$dS,nbins=10,method="content") ##create equal-sized bins (method=content)

##dS Separate into bins!! redo bin
dSs<-dS[!is.na(dS$avg_prom_meth),]
dSs$bin2<-bin(dSs$dS,nbins=10,method="content") ##create bins,equal-sized with method="content"


####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################
## PLOTTING

########################################################################################## 
## Are there differences in methylation across gene features?
ggplot(RS.fil, aes(x=region, y=pct)) +
  geom_violin(aes(fill=region), show.legend=FALSE,draw_quantiles=c(0.5)) +
  scale_x_discrete(labels=c("intergenic_nonpromoter"="Intergenic", "introns"="Introns", "exons"="Exons", "promoters"="Promoters")) +
  labs(x="", y="Percent Methylation") 

## Methylation per genome position
MTscaf<-RS.fil[RS.fil$chr=="MT",]
groups<-RS.fil[grep("group",RS.fil$chr),]
scaffs<-RS.fil[grep("scaffold_",RS.fil$chr),]

scaffs2<-scaffs
scaffs2$chr<-"scaffold"
##across chromosomes
aggs<-aggregate(pct ~ region + chr, rbind(groups,scaffs2), mean)
aggs2<-aggs
aggs2$chr<-factor(aggs2$chr,levels=c("groupI","groupII","groupIII","groupIV","groupV","groupVI","groupVII","groupVIII","groupIX","groupX","groupXI","groupXII","groupXIII","groupXIV","groupXV","groupXVI","groupXVII","groupXVIII","groupXIX","groupXX","groupXXI","scaffold"))
ggplot(aggs2, aes(x=region, y=pct,fill=region)) +
  facet_grid(~chr, switch='y') +
  geom_bar(stat = "identity") 

##########################################################################################
## CNVs 

dd<-melt(full[!is.na(full$CNVs),colnames(full) %in% c("avg_prom_meth","CNVs")])
g1<-ggboxplot(dd, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="promoter methylation",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) 

ee<-melt(full[!is.na(full$CNVs),colnames(full) %in% c("div_prom_meth","CNVs")])
ee$value<-abs(ee$value)
g3<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="promoter methylation divergence",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) 

ee<-melt(full[!is.na(full$CNVs),colnames(full) %in% c("avg_exp","CNVs")])
ee$value<-log2(ee$value)
g2<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) 

ee<-melt(full[!is.na(full$CNVs),colnames(full) %in% c("div_exp","CNVs")])
ee$value<-abs(ee$value)
g4<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="expression divergence",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) 


####################################################################################################################################################################################################################################
#########   Analyze with expression data   ######################################################################################################################################################################################
####################################################################################################################################################################################################################################

##log of 0 is Inf, so add a constant
cor.test(full_all$promoter_methylation, log2(full_all$average+0.1),na.rm=TRUE)
cor.test(full_all$genebody_methylation, log2(full_all$average+0.1),na.rm=TRUE)

pm<-ggplot(full_all, aes(x=promoter_methylation,y=log2(average))) + 
  geom_point(alpha=0.1)  + stat_smooth(method="lm") + 
  xlab("promoter methylation") + ylab("log2 expression") 
gb<-ggplot(full_all, aes(x=genebody_methylation,y=log2(average))) + 
  geom_point(alpha=0.1) + stat_smooth(method="lm") + 
  xlab("genebody methylation") + ylab("log2 expression") 

## CNVs of duplicate pairs
ee<-melt(full[,colnames(full) %in% c("avg_exp","CNVs")])
ee$value<-log2(ee$value)
a1<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) 

ee<-melt(full[,colnames(full) %in% c("avg_prom_meth","CNVs")])
ee$value<-abs(ee$value)
a2<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="promoter methylation",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ee<-melt(full[,colnames(full) %in% c("avg_gb_meth","CNVs")])
ee$value<-abs(ee$value)
a3<-ggboxplot(ee, x="CNVs", y="value",fill="CNVs",varwidth=F,notch=T,ylab="genebody methylation",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("0","1"),c("1","2"),c("0","2")),na.rm = TRUE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

### Methylation / Expression vs Categories
##By category, all genes
fmeth<-ggboxplot(full_all, x="Cat", y="promoter_methylation",fill="Cat",varwidth=F,notch=T,outlier.shape=NA,
      ylab="promoter methylation",xlab="",legend="none") +
    theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium")),na.rm = TRUE) 

fgb<-ggboxplot(full_all, x="Cat", y="genebody_methylation",fill="Cat",varwidth=F,notch=T,
      ylab="genebody methylation",xlab="",legend="none") +
  theme(axis.text.x=element_blank()) +
  stat_compare_means(comparisons = list(c("Young","Medium")),na.rm = TRUE) 

methexp2<-full_all
methexp2$average<-log2(methexp2$average)
fexp<-ggboxplot(methexp2, x="Cat", y="average",fill="Cat",varwidth=F,notch=T,
      ylab="average expression",xlab="",legend="none") +
  stat_compare_means(comparisons = list(c("Young","Medium")),na.rm = TRUE) 

##########################################################################################
##By ancestor, all genes
fmeth<-ggplot(full, aes(x=last_common_ancestor, y=abs(promoter_methylation),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  coord_cartesian(ylim = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation") + xlab("last common ancestor") 
fgb<-ggplot(methexp, aes(x=last_common_ancestor, y=(genebody_methylation),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("genebody methylation") + xlab("last common ancestor") 
fexp<-ggplot(methexp, aes(x=last_common_ancestor, y=log2(average),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("average expression") + xlab("last common ancestor") 

##########################################################################################
##Divergence
##By ancestor, all genes
fmeth2<-ggplot(full[full$prom_counts>=mincountthreshold & full$p_prom_counts>=mincountthreshold,], aes(x=last_common_ancestor, y=abs(div_prom_meth),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("promoter methylation divergence") + xlab("last common ancestor")
fgb2<-ggplot(full[full$prom_counts>=mincountthreshold & full$p_prom_counts>=mincountthreshold,], aes(x=last_common_ancestor, y=abs(div_gb_meth),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("genebody methylation divergence") + xlab("last common ancestor") 
fexp2<-ggplot(full[full$prom_counts>=mincountthreshold & full$p_prom_counts>=mincountthreshold,], aes(x=last_common_ancestor, y=abs(div_exp),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25,notch=T,show.legend=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("average expression divergence") + xlab("last common ancestor") 

##########################################################################################
##By dS, duplicates
fexp<-ggplot(dS, aes(x=bin2, y=log2(avg_exp),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  ylab("average expression") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12)) +
  xlab("dS") +
  theme(legend.position="none")
fgb<-ggplot(dS, aes(x=bin2, y=abs(genebody_methylation),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  ylab("genebody methylation") + 
  theme(legend.position="none") 
fmeth<-ggplot(dS, aes(x=bin2, y=abs(promoter_methylation),fill=bin2)) +
  geom_boxplot(notch=T,outlier.size=.25, varwidth=F) +
  ylab("promoter methylation") + 
  theme(legend.position="none") 

##################################################################################################################################################
##By ancestor, duplicate genes divergence

fexp2<-ggplot(full, aes(x=last_common_ancestor, y=abs(div_exp),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25, notch=T,show.legend=FALSE,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +
  font("ylab")+
  font("xy.text") +
  xlab("") + ylab("expression divergence") 
  
fmeth2<-ggplot(full, aes(x=last_common_ancestor, y=abs(div_prom_meth),fill=last_common_ancestor)) +
  geom_boxplot(outlier.size=.25, notch=T,show.legend=FALSE,outlier.shape=NA)
  font("ylab")+
  font("xy.text")+
  xlab("") + ylab("methylation divergence") 

##################################################################################################################################################
##expression across tissues by age

avgs<-melt(full_all[full_all$pairs=="pair",names(full_all) %in% c("Cat","last_common_ancestor","average","HK","SP","ovaries","testes","em","lm","mg","eo","X24")])
ggplot(avgs, aes(x=Cat, y=log2(value),col=Cat)) +
  geom_boxplot(outlier.size=.25, varwidth=F,notch=T) +
  labs(title="", x="", y="log2 expression") +
  theme(legend.position="top") +
  facet_grid(~variable)
