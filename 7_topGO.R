##########################################################################################
# Gene Ontology Analysis

library(topGO)

##format universe db file
biomart<-"SticklebackGO_biomart.txt" ## biomart GO annotations
biomartf<-read.table(biomart,sep="\t",header=T,stringsAsFactors=FALSE)
biomartf<-aggregate(.~Gene.stable.ID, data=biomartf, paste, collapse = ",")
db<-"Genes.IDs.db"
write.table(biomartf,file=db,row.names=F,quote=F,sep="\t",col.names=F)

##########################################################################################
##run topGO

# read custom annotation, all GO IDs separated by commas
GOmap <- readMappings(db)
# save the geneID in the universe
refset <- names(GOmap)

##input gene lists - simple list of ensembl IDs
glists<-c("YoungDuplicates.IDs","MediumDuplicates.IDs","OldDuplicates.IDs","SSDDuplicates.IDs","CNVsDuplicates.IDs","OhnoDuplicates.IDs")
for (file in glists){
  #genes in the subset, must be present in the universe
  genes<-read.table(file,header=F)
  #this extracts the 1st column that has the gene ids
  geneIDs<-genes[,1]

  #create vector describing the interesting genes with 1, the rest with 0
  genes_of_interest = factor(as.integer(refset %in% geneIDs))
  names(genes_of_interest) <- refset

  #Set the gene ontology categories to loop through
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)){
	##create the topGO data object
	tgData = new("topGOdata", ontology = ontology[i], allGenes = genes_of_interest, annot = annFUN.gene2GO, gene2GO = GOmap)
	##sets statistics 
	fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
	fisherResCor = p.adjust(score(fisherRes), method="fdr")
	weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
	weightResCor = p.adjust(score(weightRes), method="fdr")
	allRes    = GenTable(tgData, classic=fisherRes, weight=weightRes, orderBy="weight", ranksOf="classic", topNodes=100, numChar=200)
	allRes$fisher.COR = fisherResCor[allRes$GO.ID]
	allRes$weight.COR = weightResCor[allRes$GO.ID]
	write.csv(allRes, paste(file,db,"topGO.over",ontology[i],"csv",sep="."))
	}
}
##########################################################################################
# make a dot plot comparing the enriched functions between groups of genes

library(reshape2)
library(ggplot2)
library(plyr)
library(readr)

mydir="."
BPfilesG = list.files(path=mydir, pattern="*Genes.IDs.db.topGO.over.BP.csv", full.names=TRUE)
MFfilesG = list.files(path=mydir, pattern="*Genes.IDs.db.topGO.over.MF.csv", full.names=TRUE)
CCfilesG = list.files(path=mydir, pattern="*Genes.IDs.db.topGO.over.CC.csv", full.names=TRUE)

#filenames
filenames<-vapply(strsplit(BPfilesG, "\\."), `[`, 2, FUN.VALUE=character(1))###split to remove suffix
filenames<-vapply(strsplit(filenames, "\\/"), `[`, 2, FUN.VALUE=character(1))###split to remove suffix

##load all csv files together, and add names to them
csv = ldply(BPfilesG, read_csv)
csv$Genes <- rep(filenames,each=100)
csv$Cat<-"BP"
GOenrich <- csv
csv = ldply(MFfilesG, read_csv)
csv$Genes <- rep(filenames,each=100)
csv$Cat<-"MF"
GOenrich <- rbind(GOenrich,csv)
csv = ldply(CCfilesG, read_csv)
csv$Genes <- rep(filenames,each=100)
csv$Cat<-"CC"
GOenrich <- rbind(GOenrich,csv)

##only keep GO < threshold, and reorder factors
GOenrich<-GOenrich[GOenrich$weight.COR<0.05,]
GOenrich$Fold<-GOenrich$Significant/GOenrich$Expected
GOenrich<-GOenrich[order(GOenrich$Genes,-GOenrich$weight.COR),]
GOenrich$Term<-factor(GOenrich$Term, levels = rev(unique(GOenrich$Term))) ##reorder factor levels based on order in table

ggplot(GOenrich,aes(x=Genes,y=Term)) + 
  geom_point(aes(size=Fold,color=weight.COR)) +
  #scale_color_gradient2(low='lightblue',mid="lightblue", high="darkblue") +
  theme_bw() + 
  theme(axis.text.y=element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())
