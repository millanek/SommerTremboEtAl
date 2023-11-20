########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 20th Nov 2023
### Analyses about association of (log2-transformed) gene expression with the behaviour median values

########################################################################
#######   Load packages  ######################################
########################################################################


# Change the path to point to the directory containing this file
setwd("~/SommerTremboEtAl/eWAS/") 

# This is for the pgls approach
library("caper")
library("topGO"); 
# and our functions that are needed below:
source("../ExploringHAVs/GWAS_functions.R")

################################################################################
##  load data and get it into shape       ######################################
################################################################################

# Load expression data:
expression <- read.table("CichlidX_TPM_GeneExpressionMatrix_BR.txt.gz")
expressionLiver <- read.table("CichlidX_TPM_GeneExpressionMatrix_VE.txt.gz")
expressionOvary <- read.table("CichlidX_TPM_GeneExpressionMatrix_OV.txt.gz")
expressionTestis <- read.table("CichlidX_TPM_GeneExpressionMatrix_TE.txt.gz")
expressionGills <- read.table("CichlidX_TPM_GeneExpressionMatrix_GI.txt.gz")
expressionLPJ <- read.table("CichlidX_TPM_GeneExpressionMatrix_LP.txt.gz")

# the behaviour measurement
y <- read.table("../TestingForAssociation/exploratoryBehaviorMedians.txt",header=T)
# The phylogenetic tree
tanTree <- read.tree("../TestingForAssociation/TanganyikaSpeciesTree_b1.tre"); 
# Mapping between entrez gene IDs and GO terms 
onGOmap <- readMappings(file="../ExploringHAVs/onKocherGOmapping.txt") 
# genes that are near (+/-5kb) GWAS-significant SNPs and belong to significant GO terms 
uniqueSigGenes <- read.table("GWASAssociatedGenes_permutationExcluded.txt",sep="\t",header=T) 

# Preparing the expression data for further analyses 
brainExpressionOrderedLog2 <- prepareExpressionData(expression, y)
liverExpressionOrderedLog2 <- prepareExpressionData(expressionLiver, y)
GillsExpressionOrderedLog2 <- prepareExpressionData(expressionGills, y)
TestisExpressionOrderedLog2 <- prepareExpressionData(expressionTestis, y)
OvaryExpressionOrderedLog2 <- prepareExpressionData(expressionOvary, y)
LPJExpressionOrderedLog2 <- prepareExpressionData(expressionLPJ, y)

# Subset the exploratory behavior values to the species for which we have gene expression data
matchingPheno <- match(colnames(expression),y$species_abb)[-which(is.na(match(colnames(expression),y$species_abb)))]
phenoForExpression <- y[matchingPheno,]


##########################################################################################
######   looking for association of brain gene expression and exploratory behavior 
##########################################################################################

## First using standard linear model (LM)
brainExprAssoc <- doExpressionAssociation(brainExpressionOrderedLog2, phenoForExpression, tanTree)

## And some exploratory plots:
# 1) Comparing fold-change of expression against p-values for the lm method
plot(-log10(brainExprAssoc$p_GLM),abs(brainExprAssoc$foldChange),ylab="fold-change in expression (per 'unit' of behaviour)",xlab="-log10 lm eWAS p-value")
# 2) Comparing fold-change of expression against p-values for the pgls method
plot(-log10(brainExprAssoc$p_pGLS),abs(brainExprAssoc$foldChange),ylab="fold-change in expression (per 'unit' of behaviour)",xlab="-log10 pgls eWAS p-value")
# 3) Comparing p-values obtained from the lm and the pgls methods
plot(-log10(brainExprAssoc$p_GLM),-log10(brainExprAssoc$p_pGLS),ylab="-log10 pgls eWAS p-value",xlab="-log10 lm eWAS p-value")
points(-log10(brainExprAssoc$p_GLM[which(names(brainExprAssoc$p_GLM) == "100698262")]),-log10(brainExprAssoc$p_pGLS[which(names(brainExprAssoc$p_pGLS) == "100698262")]),col="red") # And where is cacng5b?


################################################################################################
######  Defining 'Interesting' genes based purely on expression and examining what they are by Gene Ontology
################################################################################################

# We said that the 'Interesting' stuff was in the top 15% of all three statistics (fold-change, lm p-val, and pgls p-val) - this gives about the same number of genes as used in the GWAS GO
ofInterestIndex <- which(-log10(ps_eWAS_LM) > quantile(-log10(ps_eWAS_LM),0.85) & abs(foldChange) > quantile(abs(foldChange),0.85) & -log10(pgls_eWASvector) > quantile(-log10(pgls_eWASvector),0.85))


# Then we do Gene Ontology enrichment analysis on these and get pretty cool results (note: the 'GeneUniverse' here are just genes expressed in the brain)
GO_LM_PGLS_expression <- doOntology.CustomAnnot(names(ps_eWAS_LM), names(ps_eWAS_LM[ofInterestIndex]), onGOmap,10)

# These are the significant terms:
sigGOsMF_eWAS <- GO_LM_PGLS_expression[[1]][as.numeric(GO_LM_PGLS_expression[[1]]$weight) < 0.05,]; sigGOIDsMF_eWAS <- sigGOsMF_eWAS$GO.ID
sigGOsCC_eWAS <- GO_LM_PGLS_expression[[2]][as.numeric(GO_LM_PGLS_expression[[2]]$weight) < 0.05,]; sigGOIDsCC_eWAS <- sigGOsCC_eWAS$GO.ID
sigGOsBP_eWAS <- GO_LM_PGLS_expression[[3]][as.numeric(GO_LM_PGLS_expression[[3]]$weight) < 0.05,]; sigGOIDsBP_eWAS <- sigGOsBP_eWAS$GO.ID

# Then we again do a 'permutation' test to confirm the validity of these results:
# 1000 permutations of GO analysis to see to what degree the results from the original analysis could appear randomly
# run_eWAS_GO_permutations(ps_eWAS_LM,ofInterestIndex, onGOmap, sigGOIDsMF_eWAS, sigGOIDsCC_eWAS, sigGOIDsBP_eWAS)   # !!!!!!! This needed to be run just once - takes hours !!!!!!!!!
# There are no terms to exclude due to permutations !!!!

# So we just output a table with all significant terms for Cytoscape:
CytoScapeMF_eWAS <- cbind(sigGOsMF_eWAS$Term, sigGOsMF_eWAS$Term, sigGOsMF_eWAS$weight, -log10(as.numeric(sigGOsMF_eWAS$weight)), rep(2,length(sigGOsMF_eWAS$Term)))
colnames(CytoScapeMF_eWAS) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScapeCC_eWAS <- cbind(sigGOsCC_eWAS$Term, sigGOsCC_eWAS$Term, sigGOsCC_eWAS$weight, -log10(as.numeric(sigGOsCC_eWAS$weight)), rep(2,length(sigGOsCC_eWAS$Term)))
colnames(CytoScapeCC_eWAS) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScapeBP_eWAS <- cbind(sigGOsBP_eWAS$Term, sigGOsBP_eWAS$Term, sigGOsBP_eWAS$weight, -log10(as.numeric(sigGOsBP_eWAS$weight)), rep(2,length(sigGOsBP_eWAS$Term)))
colnames(CytoScapeBP_eWAS) <- c("Term","Description","p.Val","FDR","Phenotype")
CytoScape_eWAS <- rbind(CytoScapeMF_eWAS, CytoScapeCC_eWAS, CytoScapeBP_eWAS)
write.table(CytoScape_eWAS, file="expression_Cytoscape_EnrichedTerms.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=T)

# Examples of how to get the genes that create the enrichment for the terms
getGenesForGO(GO_LM_PGLS_expression[[4]], names(ps_eWAS_LM[ofInterestIndex]), "GO:0005272") # sodium channel activity
getGenesForGO(GO_LM_PGLS_expression[[5]], names(ps_eWAS_LM[ofInterestIndex]), "GO:0045202") # synapse
getGenesForGO(GO_LM_PGLS_expression[[6]], names(ps_eWAS_LM[ofInterestIndex]), "GO:0007411") # axon guidance


#########################################################################################
######  Linking association of behavior to gene expression with association to genotypes  
#########################################################################################

# PART 2a) Showing that the GWAS/GO genes have on aggregate greater change in expression than all other genes

# Need to run the association for all the tissues: (Takes about 20 minutes)
# brainExprAssoc <- doExpressionAssociation(brainExpressionOrderedLog2, phenoForExpression, tanTree) (already ran above)
liverExprAssoc <- doExpressionAssociation(liverExpressionOrderedLog2, phenoForExpression, tanTree)
gillsExprAssoc <- doExpressionAssociation(GillsExpressionOrderedLog2, phenoForExpression, tanTree)
testisExprAssoc <- doExpressionAssociation(TestisExpressionOrderedLog2, phenoForExpression, tanTree)
ovaryExprAssoc <- doExpressionAssociation(OvaryExpressionOrderedLog2, phenoForExpression, tanTree)
LPJ_ExprAssoc <- doExpressionAssociation(LPJExpressionOrderedLog2, phenoForExpression, tanTree)

# Separate the fold-change vector into GWAS/GO and non-GWAS/GO genes 
GWASgenesIndex <- which(names(brainExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)
GWASgenesIndexLiver <- which(names(liverExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)
GWASgenesIndexGills <- which(names(gillsExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)
GWASgenesIndexTestis <- which(names(testisExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)
GWASgenesIndexOvary <- which(names(ovaryExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)
GWASgenesIndexLPJ <- which(names(LPJ_ExprAssoc$p_GLM) %in% uniqueSigGenes$entrezgene_id)


### And look at the difference in expression:
# First for the brain:
GWASbrain <- abs(brainExprAssoc$foldChange[GWASgenesIndex]); nonGWASbrain <- abs(brainExprAssoc$foldChange[-GWASgenesIndex])
wt.result <- wilcox.test(GWASbrain, nonGWASbrain, conf.int=T)
# and a visualization
boxplot(absFoldChangeGWASgenes, absFoldChangeOtherGenes, ylim=c(0,1.4),ylab="fold change of gene expression")
lines(c(1,2),c(1.2,1.2)); lines(c(1,1),c(1.2,1.15)); lines(c(2,2),c(1.2,1.15)); lines(c(1.5,1.5),c(1.25,1.2)); axis(1,at=c(1,2),labels=c("GWAS + GO\nassociated genes", "all other genes"),line=1,lwd=0)

# Now the other tissues:
GWASliver <- abs(liverExprAssoc$foldChange[GWASgenesIndexLiver]); nonGWASliver <- abs(liverExprAssoc$foldChange[-GWASgenesIndexLiver])
GWASgills <- abs(gillsExprAssoc$foldChange[GWASgenesIndexGills]); nonGWASgills <- abs(gillsExprAssoc$foldChange[-GWASgenesIndexGills])
GWAStestis <- abs(testisExprAssoc$foldChange[GWASgenesIndexTestis]); nonGWAStestis <- abs(testisExprAssoc$foldChange[-GWASgenesIndexTestis])
GWASovary <- abs(ovaryExprAssoc$foldChange[GWASgenesIndexOvary]); nonGWASovary <- abs(ovaryExprAssoc$foldChange[-GWASgenesIndexOvary])
GWAS_LPJ <- abs(LPJ_ExprAssoc$foldChange[GWASgenesIndexLPJ]); nonGWAS_LPJ <- abs(LPJ_ExprAssoc$foldChange[-GWASgenesIndexLPJ])

# Write out a summary table:
GeneExpressionSummmary <- makeExpressionSummaryTable(GWASbrain, nonGWASbrain,"");
GeneExpressionSummmary <- rbind(GeneExpressionSummmary, makeExpressionSummaryTable(GWASliver, nonGWASliver,"liver"))
GeneExpressionSummmary <- rbind(GeneExpressionSummmary, makeExpressionSummaryTable(GWASgills, nonGWASgills,"gills"))
GeneExpressionSummmary <- rbind(GeneExpressionSummmary, makeExpressionSummaryTable(GWAStestis, nonGWAStestis,"testis"))
GeneExpressionSummmary <- rbind(GeneExpressionSummmary, makeExpressionSummaryTable(GWASovary, nonGWASovary,"ovary"))
GeneExpressionSummmary <- rbind(GeneExpressionSummmary, makeExpressionSummaryTable(GWAS_LPJ, nonGWAS_LPJ,"LPJ"))
write.table(GeneExpressionSummmary, "GeneExpressionSummary.tsv",quote=F,sep="\t",row.names=F,col.names=T)

# Doing the Wilcox test for all the tissues to get location shift estimates and confidence intervals for each
wt.resultLiver <- wilcox.test(GWASliver, nonGWASliver, conf.int=T)
wt.resultGills <- wilcox.test(GWASgills, nonGWASgills, conf.int=T)
wt.resultTestis <- wilcox.test(GWAStestis, nonGWAStestis,, conf.int=T)
wt.resultOvary <- wilcox.test(GWASovary, nonGWASovary, conf.int=T)
wt.resultLPJ <- wilcox.test(GWAS_LPJ, nonGWAS_LPJ, conf.int=T)
shiftEstimates <- c(wt.result$estimate, wt.resultGills$estimate, wt.resultLPJ$estimate, wt.resultLiver$estimate, wt.resultOvary$estimate,  wt.resultTestis$estimate)
cisLow <- c(wt.result$conf.int[1], wt.resultGills$conf.int[1], wt.resultLPJ$conf.int[1], wt.resultLiver$conf.int[1], wt.resultOvary$conf.int[1],  wt.resultTestis$conf.int[1])
cisHigh <- c(wt.result$conf.int[2], wt.resultGills$conf.int[2], wt.resultLPJ$conf.int[2], wt.resultLiver$conf.int[2], wt.resultOvary$conf.int[2],  wt.resultTestis$conf.int[2])

par(mar=c(5.1,5.1,4.1,2.1))
plot(shiftEstimates,pch=16,ylim=c(-0.015,max(cisHigh)+0.005),ylab="difference between HAVs genes vs. other genes\nin expression fold change",xaxt='n')
for (i in 1:6) { lines(c(i,i), c(cisLow[i],cisHigh[i])) } 
abline(h=0,lty=2)
axis(1,at=1:6,labels=c("Brain", "Gills", "LPJ", "Liver", "Ovary", "Testis"))

# Another summary table as a basis for Fig. 2B
GeneExpressionMMN <- makeExpressionTableMeanMedianNumber(GWASbrain, nonGWASbrain, "brain")
GeneExpressionMMN <- rbind(GeneExpressionMMN, makeExpressionTableMeanMedianNumber(GWASgills, nonGWASgills, "gills"))
GeneExpressionMMN <- rbind(GeneExpressionMMN, makeExpressionTableMeanMedianNumber(GWAS_LPJ, nonGWAS_LPJ, "LPJ"))
GeneExpressionMMN <- rbind(GeneExpressionMMN, makeExpressionTableMeanMedianNumber(GWASliver, nonGWASliver, "liver"))
GeneExpressionMMN <- rbind(GeneExpressionMMN, makeExpressionTableMeanMedianNumber(GWAStestis, nonGWAStestis, "testis"))
GeneExpressionMMN <- rbind(GeneExpressionMMN, makeExpressionTableMeanMedianNumber(GWASovary, nonGWASovary, "ovary"))
colnames(GeneExpressionMMN) <- c("mean","median","N_genes","pc25","pc75","sd", "gene_type","tissue"); GeneExpressionMMN <- as.data.frame(GeneExpressionMMN)

plot((1:6)-0.1, GeneExpressionMMN$median[which(GeneExpressionMMN$gene_type == "HAVs")],pch=16,col="blue",ylim=c(0,0.30),ylab="fold change of gene expression",xaxt='n',xlim=c(0.8,6.2))
for (i in 1:6) { lines(c(i-0.1,i-0.1), c(GeneExpressionMMN$pc25[which(GeneExpressionMMN$gene_type == "HAVs")][i],GeneExpressionMMN$pc75[which(GeneExpressionMMN$gene_type == "HAVs")][i]), col="blue") } 
points((1:6)+0.1, GeneExpressionMMN$median[which(GeneExpressionMMN$gene_type == "other")],pch=16,col="gray50")
for (i in 1:6) { lines(c(i+0.1,i+0.1), c(GeneExpressionMMN$pc25[which(GeneExpressionMMN$gene_type == "other")][i],GeneExpressionMMN$pc75[which(GeneExpressionMMN$gene_type == "other")][i]), col="gray50") }
axis(1,at=1:6,labels=c("Brain", "Gills", "LPJ", "Liver", "Testis", "Ovary"))

