##########################################################################################
## R script for getting Genes and Paralogs based on Biomart + Ohnologs database

library(plyr)
library(data.table)
library(stringr)

setwd("Genome/")
##########################################################################################
##Gene info from Ensembl with Gene stable ID and gene biotype
geneinfo<-read.table("Gasterosteus_aculeatus.Genes.bed",header=T,sep="\t")
names(geneinfo)[1]<-"ensembl_gene_id"
##########################################################################################
##Load paralogue info through biomaRt
para <- read.table(gzfile("Gasterosteus_aculeatus.Paralogs.gz"),header=T,sep="\t",col.names = c("ensembl_gene_id", "paralogue_ensembl_gene_id", "target_query_pct_id", "query_target_pct_id", "paralogue_homology_type", "last_common_ancestor","dN","dS"))
##Average the 2 percent identity values (1 is target vs query, other is vice versa)
para$avg_pct_id <- (para$target_query_pct_id + para$query_target_pct_id)/2
##Order of LCA
para$last_common_ancestor <- factor(para$last_common_ancestor, levels=c("Gasterosteus aculeatus","Perciformes","Eupercaria","Percomorphaceae","Euacanthomorphacea", "Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii","Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta",""))
##Sort LCA by age and keep only top based on target - query pct
paraS<-para[order( para$ensembl_gene_id, para$last_common_ancestor, -(para$target_query_pct_id) ),]
paraS<-paraS[!duplicated(paraS$ensembl_gene_id),]##keeps most closely related LCA
##Keep only reciprocal top hits
for (x in 1:nrow(paraS)){
  q<-as.character(paraS$ensembl_gene_id[x])
  v<-as.character(paraS$paralogue_ensembl_gene_id[x])
  recip<-as.character(paraS[paraS$ensembl_gene_id==v,2])
  if(length(recip)){ ##if duplicate exists, without this statement there is an error
    if(q==recip){
    paraS$pairs[x]<-"pair"
    }
    else{
      paraS$pairs[x]<-"no"
    }
  }else{
    paraS$pairs[x]<-"single"
  }
}
##add biotype information and concatenate info
paraS<-merge(paraS,geneinfo)
paraS$id<-apply(cbind(as.character(paraS$ensembl_gene_id), as.character(paraS$paralogue_ensembl_gene_id)), 1, function(x) paste(sort(x), collapse=";"))

##Add Ohnologs info from http://ohnologs.curie.fr/

#intermediate q<0.01
opairs3 <- read.table(("gaculeatus.Pairs.Intermediate.3R.txt"),header=T,sep="\t")
opairs3$ohno <- "ohnolog"
opairs2 <- read.table(("gaculeatus.Pairs.Intermediate.2R.txt"),header=T,sep="\t")
opairs2$ohno <- "ohnolog"

#Remove 2R if also found in 3R
rem1<-merge(opairs2,opairs3,by.x=c("Ohno1","Ohno2"),by.y=c("Ohno1","Ohno2"))
rem2<-merge(opairs2,opairs3,by.x=c("Ohno1","Ohno2"),by.y=c("Ohno2","Ohno1"))
opairs2<-opairs2[!(opairs2$Ohno1 %in% rem1$Ohno1 & opairs2$Ohno2 %in% rem1$Ohno2),]
opairs2<-opairs2[!(opairs2$Ohno1 %in% rem2$Ohno1 & opairs2$Ohno2 %in% rem2$Ohno2),]

#Remove duplicates and merge
paraSS1<-merge(paraS,opairs3,by.x=c("ensembl_gene_id","paralogue_ensembl_gene_id","Gene.type"),by.y=c("Ohno1","Ohno2","Gene.type"))
paraSS2<-merge(paraS,opairs3,by.x=c("ensembl_gene_id","paralogue_ensembl_gene_id","Gene.type"),by.y=c("Ohno2","Ohno1","Gene.type"))
paraSS3<-merge(paraS,opairs2,by.x=c("ensembl_gene_id","paralogue_ensembl_gene_id","Gene.type"),by.y=c("Ohno1","Ohno2","Gene.type"))
paraSS4<-merge(paraS,opairs2,by.x=c("ensembl_gene_id","paralogue_ensembl_gene_id","Gene.type"),by.y=c("Ohno2","Ohno1","Gene.type"))
paraSS<-rbind(paraSS1,paraSS2,paraSS3,paraSS4)
paraSD<-merge(paraS,paraSS[,c(1:12,19)],all.x=T)
paraSD<-(unique(paraSD))

##add SSD and singletons
paraSD$ohno<-ifelse(is.na(paraSD$ohno),paraSD$pairs,paraSD$ohno)
paraSD$ohno[paraSD$ohno=="pair"]<-"SSD"
paraSD$ohno[paraSD$ohno=="no"]<-"SSD"

##Three discrete age categories + singletons
paraSD$Cat[paraSD$last_common_ancestor %in% "Gasterosteus aculeatus"]<-"Young"
paraSD$Cat[paraSD$last_common_ancestor %in% c("Perciformes","Eupercaria", "Percomorphaceae","Euacanthomorphacea","Acanthomorphata", "Euteleosteomorpha","Clupeocephala", "Osteoglossocephalai","Neopterygii","Actinopterygii")]<-"Medium"
paraSD$Cat[paraSD$last_common_ancestor %in% c("Euteleostomi","Gnathostomata","Vertebrata","Chordata","Bilateria","Opisthokonta")]<-"Old"
paraSD$Cat[paraSD$last_common_ancestor==""]<-"Singleton"
paraSD$Cat <- factor(paraSD$Cat, levels=c("Young","Medium","Old","Singleton"))

#########################################################################
##orthologs in zebrafish from Biomart
danio<-read.table("SticklebackZebrafishOrthologs.txt",header=T,sep="\t")
names(danio)<-c("ensembl_gene_id","z_ensembl","z_LCA","z_homology")
danio1<- danio[!duplicated(danio$ensembl_gene_id),]##keeps only first occurrence for many2many 
##collapse many2one
dr <- aggregate(z_ensembl ~ ensembl_gene_id, danio, paste, collapse=";") #collapse
dr2 <- aggregate(z_ensembl ~ ensembl_gene_id, danio, length) #collapse and count number of zebrafish genes
names(dr2)[2]<-"z_num"
dr<- merge(dr,dr2)
dr$z_num[dr$z_ensembl==""]<-0
dr<- merge(dr,danio1[,c(1,3,4)])# add LCA and homology

##get same information but from stickleback perspective
danio2<-read.csv("ZebrafishSticklebackOrthologs.txt",header=T,sep="\t")
names(danio2)<-c("z_ensembl","ensembl_gene_id","LCA","homology")
danio3<- danio2[!duplicated(danio2$ensembl_gene_id),]##keeps only first occurrence for many2many 
##collapse many2one
dr3 <- aggregate(ensembl_gene_id ~ z_ensembl, danio2, paste, collapse=";") #collapse
dr4 <- aggregate(ensembl_gene_id ~ z_ensembl, danio2, length) #collapse and count number of zebrafish genes
names(dr4)[2]<-"num"
dr3<- merge(dr3,dr4)
dr3$num[dr3$ensembl_gene_id==""]<-0
dr3<- merge(dr3,danio3[,c(1,3,4)])# add LCA and homology

one_on_one<-dr[dr$z_homology=="ortholog_one2one",]
many_to_many<-dr[dr$z_homology=="ortholog_many2many",]
one_to_many<-dr[dr$z_homology=="ortholog_one2many",]

##pairs in stickleback that are singletons in zebrafish
zoutgroup<-(dr3[dr3$num==2 & dr3$homology=="ortholog_one2many",])

outgroups<-danio3[danio3$z_ensembl %in% zoutgroup$z_ensembl,]
outgroups2<-danio[danio$z_homology=="ortholog_one2one",]
outgroups2<-outgroups2[,c(2,1,3,4)]
names(outgroups2)[3:4]<-c("LCA","homology")
outgroups<-rbind(outgroups,outgroups2)
write.csv(outgroups, "Stickleback_Zebrafish_one2one_one2many.csv")

##########################################################################################
##Add CNVs from Chain et al 2014
DUP<- read.table("Duplications.masked.bed.genesfull.indivs.pops.bed",header=F)
DEL<- read.table("Deletions.masked.bed.genesfull.indivs.pops.bed",header=F) 
DUP$cnv<-"dup"
DEL$cnv<-"del"
CNV<- rbind(DUP,DEL)
CNV<- CNV[,c(4,13,15,17)]
names(CNV)<- c("ensembl_gene_id","NumInd","NumPop","cnv")

##Reorder and remove repeats
CNV<- CNV[order(-CNV$NumInd),]
CNV<- CNV[!duplicated(CNV$ensembl_gene_id),]

full<-merge(paraSD,CNV,all.x=T)
full$cnv[is.na(full$cnv)]<-0
full$CNV<-ifelse(full$cnv=="del"|full$cnv=="dup",1,0)

##########################################################################################
##Add orthologs Zebrafish, and combine/sum columns
parao<-merge(full,dr)
ppairs<-merge(parao,parao,by.x=c("ensembl_gene_id","id"),by.y=c("paralogue_ensembl_gene_id","id"),all.x=T)
##sum cnvs
ppairs$CNVs<-rowSums(ppairs[names(ppairs) %in% c("CNV.x","CNV.y")])
ppairs$cnvs<-apply(ppairs[names(ppairs) %in% c("cnv.x","cnv.y")], 1, function(x) paste(sort(x), collapse=";"))
##concatenate/unique zebrafish
ppairs$zf_ensembl<-apply(cbind(as.character(ppairs$z_ensembl.x), as.character(ppairs$z_ensembl.y)), 1, function(x) paste0(sort(unique(x)), collapse=";"))
ppairs$zf_num<-str_count(ppairs$zf_ensembl, 'ENSDARG')
ppairs$zf_LCA<-apply(cbind(as.character(ppairs$z_LCA.x), as.character(ppairs$z_LCA.y)), 1, function(x) paste0(sort(unique(x)), collapse=";"))
ppairs$zf_homology<-apply(cbind(as.character(ppairs$z_homology.x), as.character(ppairs$z_homology.y)), 1, function(x) paste0(sort(unique(x)), collapse=";"))
##remove excess columns and rename 
names(ppairs)[4:14]<-unlist(tstrsplit(names(ppairs)[4:14],"\\.",keep=1))
ppairs<-ppairs[,!names(ppairs) %in% c("ensembl_gene_id.y","target_query_pct_id.y","query_target_pct_id.y","paralogue_homology_type.y","last_common_ancestor.y","dN.y","dS.y","avg_pct_id.y","pairs.y","Gene.type.y","ohno.y","Cat.y")]
write.csv(ppairs, "SticklebackGeneParalogs.v99.csv")

##remove duplicated ids and excess columns, and repeat
pairs<-ppairs[!duplicated(ppairs$id),]
pairs<-pairs[which(as.character(pairs$paralogue_ensembl_gene_id)==as.character(pairs$ensembl_gene_id.y)),]
names(pairs)[4:14]<-unlist(tstrsplit(names(pairs)[4:14],"\\.",keep=1))
pairs<-pairs[,!names(pairs) %in% c("ensembl_gene_id.y","target_query_pct_id.y","query_target_pct_id.y","paralogue_homology_type.y","last_common_ancestor.y","dN.y","dS.y","avg_pct_id.y","pairs.y","Gene.type.y","ohno.y","Cat.y")]
write.csv(pairs, "SticklebackGeneParalogPairs.v99.csv")
