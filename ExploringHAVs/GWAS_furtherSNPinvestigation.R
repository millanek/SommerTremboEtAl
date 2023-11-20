
########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 20th November 2023
### This script is for more detailed investigation of results of GWAS which were produced on the compute cluster by GWASrun.R + PGLSrun.R 
### 

########################################################################
#######   Load packages  	      ######################################
########################################################################

# Change the path to point to the directory containing this file
setwd("~/SommerTremboEtAl/ExploringHAVs")

library("biomaRt"); ensembl=useMart("ensembl",dataset="oniloticus_gene_ensembl")
library("topGO"); library("GenomicRanges"); 
# and our functions that are needed below:
source("GWAS_functions.R")


########################################################################
##  get and prepare data          ######################################
########################################################################

allPvalsGT2 <- read.table("GLM_piVals_gt0.01.txt.gz", header=T) # LM pvalues < 0.01
allPGLSPvalsGT2 <- read.table("pGLS_piVals_gt0.01.txt.gz",header=T) # PGLS pvalues < 0.01
# A file with joint LM and PGLS pvalues where both are < 0.01:  
GWAS_PGLS_GT2 <- read.table("GWAS_and_PGLSpiVals_gt0.01.txt.gz",header=T)

# Loading p-values from the simulations:
# We choose to use the Ne of 80000 which is more conservative than 20000
allSimPvalsGT2_recomb_80000 <-  read.table("simulations_Ne_80000_PGLSpiVals0.01.txt.gz",header=T)
allSimPvalsGT2_LM_recomb_80000 <-  read.table("simulations_Ne_80000_GLMpiVals0.01.txt.gz",header=T)
# create a set of SNPs which are above > 2 in both PGLS and LM for the overlap:
allSimPvalsGT2_recomb_80000 <- cbind(allSimPvalsGT2_recomb_80000, paste(allSimPvalsGT2_recomb_80000$chr, allSimPvalsGT2_recomb_80000$coord,sep="_"))
allSimPvalsGT2_LM_recomb_80000 <- cbind(allSimPvalsGT2_LM_recomb_80000, paste(allSimPvalsGT2_LM_recomb_80000$chr, allSimPvalsGT2_LM_recomb_80000$coord,sep="_"))
allSimPvalsGT2_recomb_80000_withLM <- allSimPvalsGT2_recomb_80000[which(allSimPvalsGT2_recomb_80000[,4] %in% allSimPvalsGT2_LM_recomb_80000[,4]),]
allSimPvalsGT2_LM_recomb_80000_withPGLS <- allSimPvalsGT2_LM_recomb_80000[which(allSimPvalsGT2_LM_recomb_80000[,4] %in% allSimPvalsGT2_recomb_80000_withLM[,4]),]

# Preparing the gene ontology annotation (only needed to be executed once, then saved in a file):
#GO <- getBM(c("entrezgene_id","go_id"),mart = ensembl)
#GOl <- aggregate(go_id~entrezgene_id, unique(GO), FUN=paste,sep=",")
#GOld <- cbind(GOl[[1]],GOl[[2]])
#write.table(GOld,"onKocherGOmapping.txt",sep="\t",quote=FALSE,row.names=FALSE)

# Load the gene ontology annotation, gene descriptions, genes themselves, and link all that stuff together using the NCBI "Entrez_ID" as the unique gene identifier
onGOmap <- readMappings(file="onKocherGOmapping.txt")  ## Mapping between entrez gene IDs and GO terms (got this from Ensembl, biomart, as described above)
geneDescriptions <- getBM(c("entrezgene_id","external_gene_name","description","ensembl_gene_id"),mart = ensembl)
genes <- read.table("GCF_001858045.1_ASM185804v2_singleCover.gp",header=T); genes <- genes[,c(2,3,4,5,6,7,8,9,13)]; genes$chrom <- LGtoNCBIVec_onKocher(genes$chrom)
m <- match(genes$name2, geneDescriptions$entrezgene_id); m2 <- match(genes$name, geneDescriptions$external_gene_name)  # Matching between gene descriptions and the annotation
orderedDesc <- character(0); for (i in 1:length(m)) { if (!is.na(m[i])) { orderedDesc <- c(orderedDesc, geneDescriptions$description[m[i]])  } else if (is.na(m2[i])) { orderedDesc <- c(orderedDesc, geneDescriptions$description[m[i]]) } else { orderedDesc <- c(orderedDesc, "") } }
genesWithDescription <- cbind(genes, orderedDesc)
genesOnChroms <- genes[grep("NC_*",genes$chrom,perl=T),]
# There are a few (72) genes (i.e. locations) which point to the same entrezID
# for the goseq method we need unique entrez IDs, so for conveninence we just iremove these 72 from the "gene universe"
uniqueIDgenesOnChroms <- genesOnChroms[!duplicated(genesOnChroms$name2),] 

# This one also includes the R^2 values, which we found are perfectly correlated with the p-values, so probably not needed here; if someone asks about the R2 values, they can be found here
#allPvalsGT2 <- read.table("GLM_piVals_andRsquared.txt.gz",header=T)


##########################################################################################################################
##  subsetting SNPs to look at the top 0.1% (strict) and the top 1% (permissive)       				 #####################
##########################################################################################################################

topLMsnpsPermissive <- allPvalsGT2[which(-log10(allPvalsGT2$ps) > quantile(-log10(allPvalsGT2$ps),0.99)),]
topPGLSsnpsPermissive <- allPGLSPvalsGT2[which(-log10(allPGLSPvalsGT2$ps) > quantile(-log10(allPGLSPvalsGT2$ps),0.99)),]

# And we need the overlap of PGLS and LM for the more 'permissive' ones 
topSNPsPermissive_LM_PGLS <- GWAS_PGLS_GT2[which((-log10(GWAS_PGLS_GT2$ps) > quantile(-log10(allPvalsGT2$ps),0.99)) & (-log10(GWAS_PGLS_GT2$PGLSps) > quantile(-log10(allPGLSPvalsGT2$ps),0.99))),]

# The number of variants that pass these thresholds in simulations
length(which(-log10(allSimPvalsGT2_LM_recomb_80000_withPGLS$ps) > quantile(-log10(allPvalsGT2$ps),0.99) & -log10(allSimPvalsGT2_recomb_80000_withLM$pglsGWASvector) > quantile(-log10(allPGLSPvalsGT2$ps),0.99)))


##########################################################################################################################
##  plot Fig 2a - overlaying the simulations on the LM and PGLS pVals from real data        		#####################
##########################################################################################################################

# Plots with simulated and real data:
png(paste0("Fig2a.png"),width=1000,height=1000)
plot(-log10(GWAS_PGLS_GT2$ps),-log10(GWAS_PGLS_GT2$PGLSps),xlab="GWAS (LM) p-values",ylab="PGLS p-values",cex.axis=1.5,cex.lab=1.5,pch=16,cex=0.5)
points(-log10(allSimPvalsGT2_LM_recomb_80000_withPGLS$ps),-log10(allSimPvalsGT2_recomb_80000_withLM$pglsGWASvector),col="gold1",pch=16,cex=0.5)
legend("topleft",legend=c("Real SNPs","Simulated SNPs"),col=c("black","gold1"),pch=c(1,1),cex=2)
dev.off()


##########################################################################################################################
##  investigating the 'corners' of the plot    											     		#####################
##########################################################################################################################


# What are these weird SNPs which are high in one method and not the other, in a combination of values that does not appear in the simulations:
# We will make tree-visualisations for a few selected SNPs:
# a) bottom-right corner:
GWAS_PGLS_GT2[which(-log10(GWAS_PGLS_GT2$ps) > 11 & -log10(GWAS_PGLS_GT2$PGLSps) < 2.1),]
# b) top-left-corner:
GWAS_PGLS_GT2[which(-log10(GWAS_PGLS_GT2$ps) < 2.1 & -log10(GWAS_PGLS_GT2$PGLSps) > 4),]


##########################################################################################################################
##  Finding genes near (within 5kb) the 'top' SNPs 									     			 #####################
##########################################################################################################################


# Create a GenomicRanges object (from the GenomicRanges R package) which allows us to find overlaps and relative positions between genomic regions easily
allGenesGR <- makeGRangesFromDataFrame(genesOnChroms, keep.extra.columns = TRUE,start.field="txStart", end.field="txEnd", seqnames.field="chrom")

# Combined PGLS and LM (SNPs that are high in both)
near5kbGenesPGLS_LM <- findGenesNearSNPs(topSNPsPermissive_LM_PGLS, allGenesGR, 5000)
near5kbGenesPGLS_LM_IDs <- unique(near5kbGenesPGLS_LM$name2)
# Write into a file the entrez IDs of all the genes that have a GWAS SNP (combined PGLS and LM) within 5kb of them
write.table(near5kbGenesPGLS_LM_IDs, "GWAS_allAssociatedGenes5kb_PGLSplusLM.txt",row.names=F,col.names=F,quote=F)

allGWASassociatedGenes <- getBM(c("entrezgene_id","description","external_gene_name","chromosome_name","start_position","end_position"),filters="entrezgene_id",values=near5kbGenesPGLS_LM_IDs,mart = ensembl)

allGWASassociatedGenes <- genesWithDescription[which(genesWithDescription$name2 %in% near5kbGenesPGLS_LM_IDs),]
write.table(allGWASassociatedGenes, "GWAS_allAssociatedGenes5kb_PGLSplusLM_withDetails.txt",row.names=F,col.names=F,quote=F,sep="\t")

##########################################################################################################################
##  Doing the Gene Ontology analysis with permutations								     			 #####################
##########################################################################################################################

GO_PGLS_LM_5kb <- doOntology.CustomAnnot(genesOnChroms$name2, near5kbGenesPGLS_LM_IDs, onGOmap,10)

sigGOsMF <- GO_PGLS_LM_5kb[[1]][as.numeric(GO_PGLS_LM_5kb[[1]]$weight) < 0.05,]
sigGOsCC <- GO_PGLS_LM_5kb[[2]][as.numeric(GO_PGLS_LM_5kb[[2]]$weight) < 0.05,]
sigGOsBP <- GO_PGLS_LM_5kb[[3]][as.numeric(GO_PGLS_LM_5kb[[3]]$weight) < 0.05,]

sigGOIDsMF <- sigGOsMF$GO.ID
sigGOIDsCC <- sigGOsCC$GO.ID
sigGOIDsBP <- sigGOsBP$GO.ID

# 1000 permutations of GO analysis to see to what degree the results from the original analysis could appear randomly:
#--- This needed to be done only once - the results are witten into text files at the end of this function 
doGOpermutations()

# Next, we just load the results directly from the text files to see what terms to exclude
permutationResMF <- read.table("permutationResMF_means.txt",sep="\t")
permutationResCC <- read.table("permutationResCC_means.txt",sep="\t")
permutationResBP <- read.table("permutationResBP_means.txt",sep="\t")
# Excluding terms appearing in more than 10% of the permutation runs
permutationExcludes <- rbind(permutationResMF[which(permutationResMF[,2] > 0.1),],permutationResCC[which(permutationResCC[,2] > 0.1),],permutationResBP[which(permutationResBP[,2] > 0.1),])
sigGOsMF_np <-  sigGOsMF[-which(sigGOsMF$Term %in% permutationExcludes[,1]),]
sigGOsCC_np <-  sigGOsCC[-which(sigGOsCC$Term %in% permutationExcludes[,1]),]
sigGOsBP_np <-  sigGOsBP[-which(sigGOsBP$Term %in% permutationExcludes[,1]),]

##########################################################################################################################
##  Generating the Enriched Terms table for Cytoscape								     			 #####################
##########################################################################################################################

CytoScapeMF <- cbind(sigGOsMF_np$Term, sigGOsMF_np$Term, sigGOsMF_np$weight, -log10(as.numeric(sigGOsMF_np$weight)), rep(2,length(sigGOsMF_np$Term)))	
colnames(CytoScapeMF) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScapeCC <- cbind(sigGOsCC_np$Term, sigGOsCC_np$Term, sigGOsCC_np$weight, -log10(as.numeric(sigGOsCC_np$weight)), rep(2,length(sigGOsCC_np$Term)))	
colnames(CytoScapeCC) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScapeBP <- cbind(sigGOsBP_np$Term, sigGOsBP_np$Term, sigGOsBP_np$weight, -log10(as.numeric(sigGOsBP_np$weight)), rep(2,length(sigGOsBP_np$Term)))	
colnames(CytoScapeBP) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScape <- rbind(CytoScapeMF, CytoScapeCC, CytoScapeBP)
CytoScapeTop5 <- rbind(CytoScapeMF[1:5,], CytoScapeCC[1:5,], CytoScapeBP[1:5,])
# This will need manual editing to complete the terms ending in three dots "..."
write.table(CytoScape, file="GWAS_Cytoscape_EnrichedTerms.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=T)
write.table(CytoScapeTop5, file="~/CarolinGWAS/results/GWAS_Cytoscape_EnrichedTerms_top5.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=T)
# A version of this table that has been manually edited is in /Users/milanmalinsky/CarolinGWAS/210520_EnrichedTermsForNetwork.GO_PGLS_LM_5kb_removedNonsignificantPermutationTerms.txt


##########################################################################################################################
##  Generating tables of genes associated with GO-significant terms with significant SNPs within 5kb of them	##########
##########################################################################################################################

allSigGenes <- numeric(0); i <- 1;
for (i in 1:dim(sigGOsMF_np)[1]) {
	term <- sigGOsMF_np$GO.ID[i]
	termName <- sigGOsMF_np$Term[i]

	if (i == 1) {		
		allSigGenes <- getGenesForGO(GO_PGLS_LM_5kb[[4]], near5kbGenesPGLS_LM_IDs, term)
		n <- names(allSigGenes)
		allSigGenes <- cbind(rep(termName,dim(allSigGenes)[1]), allSigGenes)
		names(allSigGenes) <- c("GO_term", n)
	} else {	
		theseSigGenes <- getGenesForGO(GO_PGLS_LM_5kb[[4]], near5kbGenesPGLS_LM_IDs, term)
		n <- names(theseSigGenes)
		theseSigGenes <- cbind(rep(termName,dim(theseSigGenes)[1]), theseSigGenes)
		names(theseSigGenes) <- c("GO_term", n)
		allSigGenes <- rbind(allSigGenes, theseSigGenes)
	}
}

for (i in 1:dim(sigGOsCC_np)[1]) {
	term <- sigGOsCC_np$GO.ID[i]
	termName <- sigGOsCC_np$Term[i]
	theseSigGenes <- getGenesForGO(GO_PGLS_LM_5kb[[5]], near5kbGenesPGLS_LM_IDs, term)
	n <- names(theseSigGenes)
	theseSigGenes <- cbind(rep(termName,dim(theseSigGenes)[1]), theseSigGenes)
	names(theseSigGenes) <- c("GO_term", n)
	allSigGenes <- rbind(allSigGenes, theseSigGenes)
}

for (i in 1:dim(sigGOsBP_np)[1]) {
	term <- sigGOsBP_np$GO.ID[i]
	termName <- sigGOsBP_np$Term[i]
	
	theseSigGenes <- getGenesForGO(GO_PGLS_LM_5kb[[6]], near5kbGenesPGLS_LM_IDs, term)
	n <- names(theseSigGenes)
	theseSigGenes <- cbind(rep(termName,dim(theseSigGenes)[1]), theseSigGenes)
	names(theseSigGenes) <- c("GO_term", n)
	allSigGenes <- rbind(allSigGenes, theseSigGenes)
}

write.table(allSigGenes,"GWAS_allSigGenes_withTerms_permutationExcluded.txt",sep="\t",quote=F,row.names=F)

uniqueSigGenes <- unique(allSigGenes[,-1])
write.table(uniqueSigGenes,"GWAS_uniqueSigGenes_permutationExcluded.txt",sep="\t",quote=F,row.names=F)


##############################################################################################
######   Making Fig 2D			 						 			 #####################
##############################################################################################


load("GLM_results_NC_031969.gwas",verbose=T) # The GLM objects are loaded for a particular chromosome

# Get the lm result for a specific SNP and assign it to the 'a' variable
a <- GWASvector.noMiss[[which(names(GWASvector.noMiss) == "1472288")]] # Getting the distribution of alleles across species for that SNP
y <- read.table("../TestingForAssociation/exploratoryBehaviorMedians.txt",header=T)
y2 <- y[-1,]; y3 <- y2[order(y2$species_abb),]
SNP_species_info <- cbind(a$model,y3$species_abb); names(SNP_species_info) <- c("expl","gen_highest","taxa")

library(vioplot) # You may need to do install.packages("vioplot") first
pdf("Fig_2D.pdf", width=7, height=7)
vioplot(y3$median_exploration[which(a$model[,2] == 1)], y3$median_exploration[which(a$model[,2] == 0.5)], y3$median_exploration[which(a$model[,2] == -1)],names=c("genotypes: 1/1 1/1", "genotypes: 1/1 1/0", "genotypes: 0/0 0/0"),ylab="exploration")
dev.off()
