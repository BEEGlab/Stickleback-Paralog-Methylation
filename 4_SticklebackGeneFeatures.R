##########################################################################################
#R script to determine promoter-wide methylation of each gene and gene-wide methylation of each gene for an individual
#Assumes already downloaded/formatted gene lists from Ensembl
#usage: Rscript 4_SticklebackGeneFeatures.R InputFile.bismark.cov.gz

library(data.table) ##allows tstrsplit

bedfile="Gasterosteus_aculeatus.Genes.bed" #biomaRt file
bedfeature="Gasterosteus_aculeatus.gff3.genePred.bed" ###formatted GFF, from Ensembl, see CGproportionsFromFASTA.R script
CGfile="CG.stickleback.bed"

### Load annotations ... ASSIGN ENSEMBL TRANSCRIPT AND GENE IDs
# the feature annotation step adds gene IDs for sites in gene features, but this step should get all the intergenic nonpromoter sites in addition to genic sites

bed<-read.table(bedfile,header=T,sep="\t")
colnames(bed)[c(1:5,9)] <- c("seqname", "start", "end", "ensembl_gene_id","ensembl_transcript_id","biotype")
bed$start<-bed$start-1 ### make 0-based

###keep only 1 transcript per gene, the earlierst occurrence
bed$earliest<-ifelse(bed$Strand=="1",bed$Transcript.start..bp.,bed$Strand*bed$Transcript.end..bp.)
bed<-bed[order(bed$ensembl_gene_id,(bed$earliest)),] ##sort by earliest, then only keep first occurence
bed<-bed[!duplicated(bed$ensembl_gene_id),] ##gets rid of those with multiple transcripts
bed <- bed[,c(1:5,9)] ##by earliest
bed$length<-bed$end-bed$start ##gene length
write.csv(bed, paste(bedfile,"csv",sep=".")) ##gene information 

##########################################################################################
##Only  keep one transcript per gene from the GFF file
features<-read.table(bedfeature,sep="\t")
#rename transcript
features$V4<-unlist(tstrsplit(as.character(features$V4),"transcript:",keep=2))
#only keep lines with transcripts of interest (one per gene)
features<-features[features$V4 %in% bed$ensembl_transcript_id,]
write.table(features, paste(bedfeature,"filtered",sep="."),row.names=F,col.names=F,quote=F,sep="\t") ## features
