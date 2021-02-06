# SticklebackParalogMethylation
This repository contains R scripts used to analyze duplicate gene DNA methylation in threespine sticklebacks. R v3.6 was used.

## Files and main output files

### 1. Paralogs_Biomart_Ohnologs.R
Get data from biomart and determine orthology / paralogy for each gene.   
  * write.csv(outgroups, "Stickleback_Zebrafish_one2one_one2many.csv") ##list of orthologs with zebrafish.  
  * write.csv(ppairs, "SticklebackGeneParalogs.v99.csv") ##list of paralogs. 
  * write.csv(pairs, "SticklebackGeneParalogPairs.v99.csv") ##list of paralog pairs. 

### 2. ExpressionAnalysis.R
Get expression data and add to paralog/ortholog info. 
  * write.csv(sti, "Transcriptome.TPM.Stickleback.csv") ###Expression in TPM of several stickleback tissues /dev. 
  * write.csv(zebra, "Transcriptome.TPM.Zebrafish.csv") ###Expression in TPM of several zebrafish tissues /dev. 
  * write.csv(e, "Transcriptome.TPM.dup.div.exp.csv") ###Info on stickleback gene expression and duplicate gene expression divergence. 

### 3. CGproportionsFromFASTA.R
Get CpG site distribution and CpG methylation across the genome, and calculate density around TSS. 

### 4. SticklebackGeneFeatures.R
Get gene features into files for processing. 
  * write.csv(bed, paste(bedfile,"csv",sep=".")) ## gene information. 
  * write.table(features, paste(bedfeature,"filtered",sep="."),row.names=F,col.names=F,quote=F,sep="\t") ## features. 

### 5. Methylation_Individual_Analysis_Revised.R
Calculate methylation stats and trends per individual. 

### 6. Methylation_CompareIndividuals_Revised.R
Calculate differences among individuals and averaged across individuals. 
  * write.csv(allg, "Combined.gene.methylation.csv",row.names=F) ##combined methylation of each gene from each indiv. 
  * write.csv(dupg, "Combined.dup.methylation.csv",row.names=F) ##combined methylation of each paralog from each indiv. 
  * write.csv(a, "Concatenated.gene.methylation.csv",row.names=F) ##concatenated information for all genes. 
  * write.csv(d, "Concatenated.dup.methylation.csv",row.names=F) ##concatenated information for duplicate genes. 
  * write.csv(ap, "Concatenated.gene.methylation.indivs.csv",row.names=F) ##concatenated information for all genes per indiv. 
  * write.csv(dp, "Concatenated.dup.methylation.indivs.csv",row.names=F) ##concatenated information for duplicate genes per indiv. 
  * write.csv(Amethf, "Merged.methylation.filtered.all2.csv",row.names=F) ##combined information for genes genes across all samples. 
  * write.csv(Dmethf, "Merged.methylation.filtered.dups2.csv",row.names=F) ##combined information for duplicate genes across all samples. 

### 7. topGO.R
Gene Ontology Enrichment analyses

