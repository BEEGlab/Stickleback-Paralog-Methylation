##########################################################################################
##R script to calculate CG proportions in a FASTA file and methylation around TSS

library(Biostrings)
library(seqinr)
library(stringr)

##directory containing FASTA sequences of each contig/scaffold in the genome
file = list.files(path = "FASTA")
##loop through files to calculate CpG sites
for (k in 1:length(file)){
 y<- sprintf("FASTA/%s",file[k])
 readdata=readDNAStringSet(y)  
 seq = paste(readdata)
 nam = strsplit(names(readdata)," ")[[1]][1]
 val = str_locate_all(pattern ='CG', seq)
 if(length(val[[1]])>0){
  final <- data.frame(nam, val, "CG")
  write.table(final,"CG.stickleback.bed",quote=F, row=F, append = TRUE, col.names = F)
 }
}

####################################################################################################################################
##Determine CG and methylation around Transcription Start Sites

library(GenomicRanges)
library(genomation)
library(data.table)

##input files
CG<-read.table("CG.stickleback.bed") ### CG File created above from FASTA file
methfile="Combined.bismark.cov.gz" ###methylation file after processing in bismark
## NEED TO CONVERT ENSEMBL DATA TO UCSC-style BED, using the program gff3ToGenePred and modifying transcript names
#/usr/local/bin/gff3ToGenePred Gasterosteus_aculeatus.gff3 Gasterosteus_aculeatus.gff3.genePred
#/usr/local/bin/genePredToBed Gasterosteus_aculeatus.gff3.genePred Gasterosteus_aculeatus.gff3.genePred.bed
gfffile="Gasterosteus_aculeatus.gff3.genePred.bed" ###GFF, from Ensembl

##filter parameters
depth_per_site<- 10
quantile_for_filter <- .999 

####################################################################################################################################
# Get features
feat=readTranscriptFeatures(gfffile,
                              remove.unusual = FALSE,unique.prom=FALSE,
                              up.flank = 3000, down.flank = 3000)
prom=feat$promoters # get promoters from the features
feat@unlistData@elementMetadata@listData$name<-unlist(tstrsplit(feat@unlistData@elementMetadata@listData$name,"transcript:",keep=2)) ##rename 

####################################################################################################################################
# Get CG
names(CG)<-c("chr","start","end","CG")
CG.gr <- as(CG, "GRanges") # GRanges object needed for findOverlap()

# Get CG content around TSS
sm=ScoreMatrix(CG.gr,prom,strand.aware = TRUE,is.noCovNA=TRUE)
pdf("CG.TSS.pdf")
plotMeta(sm, profile.names = "CG content", xcoords = c(-3000,3000),
         ylab="CpG",dispersion = "se",
         xlab="bases around TSS")
dev.off()

####################################################################################################################################
# Get methylation values around TSSes, needs GRanges
original <- read.table(gzfile(paste("..",methfile,sep="/")),sep="\t",header=F) ## ENSEMBL
names(original)<-c("chr","start","end","pct","num","sites")

# create coverage column and filter by coverage
original <- original[original$sites <= quantile(original$sites, quantile_for_filter), ] # making sure that taking out the lower ones does not interfere w/ 99.9th percentile calculation here
RS.gr <- as(original, "GRanges") # GRanges object needed for findOverlap()

sm=ScoreMatrix(RS.gr,prom,strand.aware = TRUE, weight.col="pct",is.noCovNA=TRUE)

# Get methyalation around TSS
pdf(paste(methfile,"Methyation.TSS.pdf",sep="."))
plotMeta(sm, profile.names = "Methylation", xcoords = c(-3000,3000),
         ylab="Methylation",dispersion = "se",
         xlab="bases around TSS")
dev.off()

####################################################################################################################################
