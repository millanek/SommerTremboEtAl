########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 20th Nov 2023
### Exploring the gene expression of cacng genes across tissues

################################################################################
##  load data and get it into shape       ######################################
################################################################################

# Change the path to point to the directory containing this file
setwd("~/SommerTremboEtAl/cacngExploration/") 

# Load expression data:
expression <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_BR.txt.gz")
expressionLiver <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_VE.txt.gz")
expressionOvary <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_OV.txt.gz")
expressionTestis <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_TE.txt.gz")
expressionGills <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_GI.txt.gz")
expressionLPJ <- read.table("../eWAS/CichlidX_TPM_GeneExpressionMatrix_LP.txt.gz")

# Preparing the expression data for further analyses 
brainExpressionOrderedLog2 <- prepareExpressionData(expression, y, exclZeros = F)
liverExpressionOrderedLog2 <- prepareExpressionData(expressionLiver, y, exclZeros = F)
GillsExpressionOrderedLog2 <- prepareExpressionData(expressionGills, y, exclZeros = F)
TestisExpressionOrderedLog2 <- prepareExpressionData(expressionTestis, y, exclZeros = F)
OvaryExpressionOrderedLog2 <- prepareExpressionData(expressionOvary, y, exclZeros = F)
LPJExpressionOrderedLog2 <- prepareExpressionData(expressionLPJ, y, exclZeros = F)

# the behaviour measurement
y <- read.table("../TestingForAssociation/exploratoryBehaviorMedians.txt",header=T)
y2 <- y[-1,]; y3 <- y2[order(y2$species_abb),]
hiSpecies <- y3$species_abb[which(y3$median_exploration > quantile(y3$median_exploration,0.50))]
loSpecies <- y3$species_abb[which(y3$median_exploration <= quantile(y3$median_exploration,0.50))]

################################################################################
##  look at cacng gene expression         ######################################
################################################################################

meanExpressionValsBrain <- getMeanExpressionWithSubgroups(brainExpressionOrderedLog2, hiSpecies, loSpecies)
meanExpressionValsLiver <- getMeanExpressionWithSubgroups(liverExpressionOrderedLog2, hiSpecies, loSpecies)
meanExpressionValsGills <- getMeanExpressionWithSubgroups(GillsExpressionOrderedLog2, hiSpecies, loSpecies)
meanExpressionValsTestis <- getMeanExpressionWithSubgroups(TestisExpressionOrderedLog2, hiSpecies, loSpecies)
meanExpressionValsOvary <- getMeanExpressionWithSubgroups(OvaryExpressionOrderedLog2, hiSpecies, loSpecies)
meanExpressionValsLPJ <- getMeanExpressionWithSubgroups(LPJExpressionOrderedLog2, hiSpecies, loSpecies)

# separating species according to the cacng5b SNP
SNP_species_info <- read.table("~/CarolinGWAS/pointOfBeauty_LG4_1472288.txt",header=TRUE)
cacng5bHiSpecies <- SNP_species_info$taxa[which(SNP_species_info$gen_highest == -1)]
cacng5bLoSpecies <- SNP_species_info$taxa[which(SNP_species_info$gen_highest == 1)]
# this can be used as a switch to either use exploration mesurements or the gacng5b genotype for separation of species for the cacng barplots in part3
hiSpecies <- cacng5bHiSpecies; loSpecies <- cacng5bLoSpecies

# entrezIDs of cacng genes:
cacng_genes <- read.table("cacngEntrezIDs.txt",header=F,sep="\t"); names(cacng_genes) <- c("gene","entrezID")
cacng_gene_positions <- findCacngs(cacng_genes, expression);
cacng_genes <- cbind(cacng_genes, cacng_gene_positions)

# Fig. S11
tissues <- c("brain","liver","ovary","testis","gills","LPJs")
par(mfrow=c(3,2))
barplot(meanExpressionValsBrain$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[1],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(meanExpressionValsLiver$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[2],ylim=c(0,7))
barplot(meanExpressionValsOvary$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[3],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(meanExpressionValsTestis$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[4],ylim=c(0,7))
barplot(meanExpressionValsGills$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[5],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(meanExpressionValsLPJ$mean[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[6],ylim=c(0,7))


# Also separated by high / low exploratory behavior (not very informative; not in the paper)
par(mfrow=c(3,2))
barplot(t(cbind(meanExpressionValsBrain$meanGroup1[cacng_gene_positions],meanExpressionValsBrain$meanGroup2[cacng_gene_positions])),beside=TRUE, names= cacng_genes$gene,las=2,main= tissues[1],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(meanExpressionValsLiver$meanGroup1[cacng_gene_positions],meanExpressionValsLiver$meanGroup2[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[2],ylim=c(0,7))
legend("topright",c("high exploration species","low exploration species"),col=c("grey10","grey70"),pch=15,cex=1.0)
barplot(t(cbind(meanExpressionValsOvary$meanGroup1[cacng_gene_positions],meanExpressionValsOvary$meanGroup2[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[3],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(meanExpressionValsTestis$meanGroup1[cacng_gene_positions],meanExpressionValsTestis$meanGroup2[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[4],ylim=c(0,7))
barplot(t(cbind(meanExpressionValsGills$meanGroup1[cacng_gene_positions],meanExpressionValsGills$meanGroup2[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[5],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(meanExpressionValsLPJ$meanGroup1[cacng_gene_positions],meanExpressionValsLPJ$meanGroup2[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[6],ylim=c(0,7))
