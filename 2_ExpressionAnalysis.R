##########################################################################################
## R script for expression analysis. Input expression datasets from featurecounts

library(ggplot2)
library(plyr)
library(reshape)

################# TPM  ####################
### Make TPM from read counts
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

genelength<-read.table("GeneLengths.txt",header=TRUE) ##list of stickleback gene lengths for TPM estimation

##########################################################################################
##Developmental transcriptomes from Kaitetzidou et al
setwd("Kaitetzidou_featureCounts/")
devinfo<-read.csv("Kaitetzidou_2019_SraRunTable.txt")
devinfo<-devinfo[,names(devinfo) %in% c("Run","dev_stage")]
devinfo$cat<-c("5","5","4","4","5","1","1","2","1","2","2","3","3","4","3") ##assign sample category

## read in featurecounts data using loop
files <- list.files(pattern=".tabular")
dev <- NULL
for (f in files) {
  dat <- read.csv(f, header=T, sep="\t", na.strings="")
  if(is.null(dev)){
    dev <-dat
  }else{
    dev <- merge(dev, dat, by="Geneid", all=TRUE)
  } 
}

##Get means
devm<-as.data.frame(dev$Geneid)
devm<-cbind(devm,"em"=rowMeans((dev[c(7,8,10)])))
devm<-cbind(devm,"lm"=rowMeans(dev[c(9,11,12)]))
devm<-cbind(devm,"mg"=rowMeans(dev[c(13,14,16)]))
devm<-cbind(devm,"eo"=rowMeans(dev[c(4,5,15)]))
devm<-cbind(devm,"24"=rowMeans(dev[c(2,3,6)]))
names(devm)[1]<-"Geneid"

##Calculate TPM
devm<-merge(devm,genelength)
tpms_ga_dev <- as.data.frame(apply(devm[,2:(ncol(devm)-1)], 2, function(x) tpm(x, devm$Length)))
row.names(tpms_ga_dev)<-as.character(devm[,1])

##########################################################################################
##Gonad data from Schmitz et al
setwd("Schmitz_featureCounts/")
goninfo<-read.csv("Gonads_SraRunTable.txt")
goninfo<-goninfo[,names(goninfo) %in% c("Run","tissue_type")]

## read data using loop
files <- list.files(pattern=".tabular")
gon <- NULL
for (f in files) {
  dat <- read.csv(f, header=T, sep="\t", na.strings="")
  if(is.null(gon)){
    gon <-dat
  }else{
    gon <- merge(gon, dat, by="Geneid", all=TRUE)
  } 
}

###Get means
ov<-rowMeans(gon[names(gon) %in% as.character(goninfo$Run[goninfo$tissue_type=="Ovaries"])])
te<-rowMeans(gon[names(gon) %in% as.character(goninfo$Run[goninfo$tissue_type=="Testes"])])
gonm<-as.data.frame(gon$Geneid)
gonm<-cbind(gonm,ovaries=as.numeric(ov),testes=as.numeric(te))
names(gonm)[1]<-"Geneid"

##Calculate TPM
gonm<-merge(gonm,genelength)
tpms_ga_gon <- as.data.frame(apply(gonm[,2:(ncol(gonm)-1)], 2, function(x) tpm(x, gonm$Length)))
row.names(tpms_ga_gon)<-as.character(gonm[,1])

##################################################################################################
##Head Kidney and Spleen data from Huang et al
setwd("Huang_featureCounts/")
hksinfo<-read.csv("SpleenHK_SraRunTable.txt")
hksinfo<-hksinfo[,names(hksinfo) %in% c("Run","title","tissue_type","Alias")]
###Replace levels/factor with short name
hksinfo$tissue<-revalue(hksinfo$tissue_type,(c("Head kidney tissue" = "HK", "Spleen tissue" = "SP")))
hksinfo$name<-paste(hksinfo$tissue,hksinfo$Alias)

## read data using loop
files <- list.files(pattern=".tabular")
hks <- NULL
for (f in files) {
  dat <- read.csv(f, header=T, sep="\t", na.strings="")
  if(is.null(hks)){
    hks <-dat
  }else{
    hks <- merge(hks, dat, by="Geneid", all=TRUE)
  } 
}

###Add together results from replicates, loop through each sample name
hkss <- NULL
for (f in sort(unique(hksinfo$name))){
  dat<-rowSums(hks[names(hks) %in% as.character(hksinfo$Run[hksinfo$name==f])])
  if(is.null(hkss)){
    hkss <- dat
  }else{
    hkss <- cbind(hkss, dat)
  } 
}
hkss<-as.data.frame(hkss)
hkss<-cbind(hks$Geneid,hkss)
names(hkss)<-c("Geneid",sort(unique(hksinfo$name)))

###Calculate means by tissue
hk<-rowMeans(hkss[names(hkss) %in% as.character(hksinfo$name[hksinfo$tissue=="HK"])])
sp<-rowMeans(hkss[names(hkss) %in% as.character(hksinfo$name[hksinfo$tissue=="SP"])])

hksm<-as.data.frame(hkss$Geneid)
hksm<-cbind(hksm,HK=as.numeric(hk),SP=as.numeric(sp))
names(hksm)[1]<-"Geneid"

##Calculate TPM
hksm<-merge(hksm,genelength)
tpms_ga_hks <- as.data.frame(apply(hksm[,2:(ncol(hksm)-1)], 2, function(x) tpm(x, hksm$Length)))
row.names(tpms_ga_hks)<-as.character(hksm[,1])

##################################################################################################
##################################################################################################
##Merge stickleback data, TPM
sti<-merge(tpms_ga_hks,tpms_ga_gon,by="row.names")
row.names(sti)<-sti[,1]
sti<-sti[,-1]
sti<-merge(sti,tpms_ga_dev,by="row.names")
names(sti)[1]<-"Geneid"


##################################################################################################
#######Zebrafish data
##from Gene Expression Atlas project, already in TPM
setwd("zebrafish")
zf<-read.table("Zebrafish_E-ERAD-475/E-ERAD-475-tpms.tsv",sep="\t",header=T) 
## get mean of replicates
zfm<-cbind(zf[,1:2],apply(zf[,3:ncol(zf)],2,function (z) sapply(strsplit(as.character(z), ",", fixed=T), function(x) mean(as.numeric(x)))))
names(zfm)[1]<-"Geneid"

##################################################################################################
##Zebrafish adult tissue data from Pasquier et al
setwd("SRP044781/featureCounts/")
zebinfo<-read.csv("Zfish_SraRunTable.txt")
zebinfo<-zebinfo[,names(zebinfo) %in% c("Run","tissue")]
zlen<-read.table("SRR1524241.tabular",header=TRUE)
  
files <- list.files(pattern=".tabular")
## read data using loop
zeb <- NULL
for (f in files) {
  dat <- read.csv(f, header=T, sep="\t", na.strings="")
  if(is.null(zeb)){
    zeb <-dat
  }else{
    zeb <- merge(zeb, dat, by="Geneid", all=TRUE)
  } 
}

###Rename columns
for (f in 2:ncol(zeb)){
  n<-names(zeb[f])
  names(zeb)[f]<-as.character(zebinfo$tissue[zebinfo$Run==n])
}

##Calculate TPM
zebm<-merge(zeb,zlen)
tpms_zebm <- as.data.frame(apply(zebm[,2:(ncol(zebm)-1)], 2, function(x) tpm(x, zebm$Length)))
row.names(tpms_zebm)<-as.character(zebm[,1])

##################################################################################################
##################################################################################################
#Merge zebrafish data, TPM
zebra<-merge(zfm,tpms_zebm,by.x="Geneid",by.y="row.names")

##################################################################################################
##################################################################################################
########################   Calculate expression characteristics     ##############################
##################################################################################################
##################################################################################################

#Here we use t() to transpose, and then cast the result as a data frame. Include any summaries in here as you wish
summarizeExp <- function(eframe){
  as.data.frame(t(apply(eframe, 1, function(v){
    average=mean(v) ##mean expression 
  })))
}

### Add summary columns, in this case mean expression across all samples
sti[,2:(ncol(sti)+1)]<-summarizeExp(sti[,2:ncol(sti)]) 
zebra[,3:(ncol(zebra)+1)]<-summarizeExp(zebra[,3:ncol(zebra)])
write.csv(sti, "Transcriptome.TPM.Stickleback.csv")
write.csv(zebra, "Transcriptome.TPM.Zebrafish.csv")

##########################################################################################
##########################################################################################
##########################################################################################
## Orthologs stickleback and zebra

setwd("Genome/")
orthopara<-read.csv("SticklebackGeneParalogs.v99.csv",row.names = 1)
## LCA -Few Eupercaria and Perciformes (<15), group with Percomorphaceae, re-configure factor
orthopara$last_common_ancestor[orthopara$last_common_ancestor=="Eupercaria" | orthopara$last_common_ancestor=="Perciformes"] <- "Percomorphaceae"
levels(orthopara$last_common_ancestor) <- c(levels(orthopara$last_common_ancestor),"Singleton")
orthopara$last_common_ancestor[orthopara$last_common_ancestor==""] <- "Singleton"
orthopara$last_common_ancestor <- factor(orthopara$last_common_ancestor, levels=c("Gasterosteus aculeatus","Percomorphaceae","Euacanthomorphacea", "Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta","Singleton"))

##########################################################################################
# expression divergence
ex<-merge(orthopara,sti,by.x="ensembl_gene_id",by.y="Geneid",all.x=T)
###Now add paralog on same row
ex2<-sti
colnames(ex2)<-paste("p",colnames(ex2),sep="_")
ex<-merge(ex,ex2,by.x="paralogue_ensembl_gene_id",by.y="p_Geneid",all.x=T)
e<-merge(orthopara,ex) ##just return to original order by Geneid

#add average and divergence expression of duplicates
e$avg_exp <- (e$average + e$p_average) / 2 
e$div_exp <- (e$average - e$p_average) / (e$average + e$p_average)
write.csv(e, "Transcriptome.TPM.dup.div.exp.csv")
