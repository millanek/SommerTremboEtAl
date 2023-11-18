

############################################################
######   PART 3 - looking specifically at cacng genes  
############################################################

# Load expression data:
expression <- read.table("CichlidX_TPM_GeneExpressionMatrix_BR.txt.gz")
expressionLiver <- read.table("CichlidX_TPM_GeneExpressionMatrix_VE.txt.gz")
expressionOvary <- read.table("CichlidX_TPM_GeneExpressionMatrix_OV.txt.gz")
expressionTestis <- read.table("CichlidX_TPM_GeneExpressionMatrix_TE.txt.gz")
expressionGills <- read.table("CichlidX_TPM_GeneExpressionMatrix_GI.txt.gz")
expressionLPJ <- read.table("CichlidX_TPM_GeneExpressionMatrix_LP.txt.gz")



matchingExpression <- match(y3$species_abb,colnames(expression))[-which(is.na(match(y3$species_abb,colnames(expression))))]
matchingExpressionBrainHi <- match(hiSpecies,colnames(expression))[-which(is.na(match(hiSpecies,colnames(expression))))]
matchingExpressionBrainLo <- match(loSpecies,colnames(expression))[-which(is.na(match(loSpecies,colnames(expression))))]
brainExpression <- expression[,matchingExpression]
expressionBrainHi <- expression[, matchingExpressionBrainHi]; expressionBrainLo <- expression[, matchingExpressionBrainLo]


tissues <- c("brain","liver","ovary","testis","gills","LPJs")
cacng_gene_positions <- numeric(0)
for (i in 1:dim(cacng_genes)[1]) {
cacng_gene_positions <- c(cacng_gene_positions,which(rownames(expression) == cacng_genes$entrezID[i]))
} 
cacng_genes <- cbind(cacng_genes, cacng_gene_positions)

100700027

expressionLiver <- read.table("~/CarolinGWAS/inputs/CichlidX_TPM_GeneExpressionMatrix_VE.txt"); 
matchingExpressionLiver <- match(y3$species_abb,colnames(expressionLiver))[-which(is.na(match(y3$species_abb,colnames(expressionLiver))))]
matchingExpressionLiverHi <- match(hiSpecies,colnames(expressionLiver))[-which(is.na(match(hiSpecies,colnames(expressionLiver))))]
matchingExpressionLiverLo <- match(loSpecies,colnames(expressionLiver))[-which(is.na(match(loSpecies,colnames(expressionLiver))))]
expressionLiverHi <- expressionLiver[, matchingExpressionLiverHi]; expressionLiverLo <- expressionLiver[, matchingExpressionLiverLo]
expressionLiver <- expressionLiver[, matchingExpressionLiver]
expressionOvary <- read.table("~/CarolinGWAS/inputs/CichlidX_TPM_GeneExpressionMatrix_OV.txt"); 
matchingExpressionOvary <- match(y3$species_abb,colnames(expressionOvary))[-which(is.na(match(y3$species_abb,colnames(expressionOvary))))]
matchingExpressionOvaryHi <- match(hiSpecies,colnames(expressionOvary))[-which(is.na(match(hiSpecies,colnames(expressionOvary))))]
matchingExpressionOvaryLo <- match(loSpecies,colnames(expressionOvary))[-which(is.na(match(loSpecies,colnames(expressionOvary))))]
expressionOvaryHi <- expressionOvary[, matchingExpressionOvaryHi]; expressionOvaryLo <- expressionOvary[, matchingExpressionOvaryLo]
expressionOvary <- expressionOvary[,matchingExpression]
expressionTestis <- read.table("~/CarolinGWAS/inputs/CichlidX_TPM_GeneExpressionMatrix_TE.txt"); 
matchingExpressionTestis <- match(y3$species_abb,colnames(expressionTestis))[-which(is.na(match(y3$species_abb,colnames(expressionTestis))))]
matchingExpressionTestisHi <- match(hiSpecies,colnames(expressionTestis))[-which(is.na(match(hiSpecies,colnames(expressionTestis))))]
matchingExpressionTestisLo <- match(loSpecies,colnames(expressionTestis))[-which(is.na(match(loSpecies,colnames(expressionTestis))))]
expressionTestisHi <- expressionTestis[, matchingExpressionTestisHi]; expressionTestisLo <- expressionTestis[, matchingExpressionTestisLo]
expressionTestis <- expressionTestis[,matchingExpression]
expressionGills <- read.table("~/CarolinGWAS/inputs/CichlidX_TPM_GeneExpressionMatrix_GI.txt"); 
matchingExpressionGills <- match(y3$species_abb,colnames(expressionGills))[-which(is.na(match(y3$species_abb,colnames(expressionGills))))]
matchingExpressionGillsHi <- match(hiSpecies,colnames(expressionGills))[-which(is.na(match(hiSpecies,colnames(expressionGills))))]
matchingExpressionGillsLo <- match(loSpecies,colnames(expressionGills))[-which(is.na(match(loSpecies,colnames(expressionGills))))]
expressionGillsHi <- expressionGills[, matchingExpressionGillsHi]; expressionGillsLo <- expressionGills[, matchingExpressionGillsLo]
expressionGills <- expressionGills[,matchingExpression]
expressionLPJ <- read.table("~/CarolinGWAS/inputs/CichlidX_TPM_GeneExpressionMatrix_LP.txt"); 
matchingExpressionLPJ <- match(y3$species_abb,colnames(expressionLPJ))[-which(is.na(match(y3$species_abb,colnames(expressionLPJ))))]
matchingExpressionLPJHi <- match(hiSpecies,colnames(expressionLPJ))[-which(is.na(match(hiSpecies,colnames(expressionLPJ))))]
matchingExpressionLPJLo <- match(loSpecies,colnames(expressionLPJ))[-which(is.na(match(loSpecies,colnames(expressionLPJ))))]
expressionLPJHi <- expressionLPJ[, matchingExpressionLPJHi]; expressionLPJLo <- expressionLPJ[, matchingExpressionLPJLo]
expressionLPJ <- expressionLPJ[, matchingExpressionLPJ]


expLogAndMeans <- function(exp) {
	expLog <- log2(exp +1)
	expression_means <- apply(expLog,1,mean)
	return(expression_means)
}

expression_means <- expLogAndMeans(brainExpression); expression_meansHi <- expLogAndMeans(expressionBrainHi); expression_meansLo <- expLogAndMeans(expressionBrainLo)
expressionLiver_means <- expLogAndMeans(expressionLiver); expressionLiver_meansHi <- expLogAndMeans(expressionLiverHi); expressionLiver_meansLo <- expLogAndMeans(expressionLiverLo)
expressionOvary_means <- expLogAndMeans(expressionOvary); expressionOvary_meansHi <- expLogAndMeans(expressionOvaryHi); expressionOvary_meansLo <- expLogAndMeans(expressionOvaryLo)
expressionTestis_means <- expLogAndMeans(expressionTestis); expressionTestis_meansHi <- expLogAndMeans(expressionTestisHi); expressionTestis_meansLo <- expLogAndMeans(expressionTestisLo)
expressionGills_means <- expLogAndMeans(expressionGills); expressionGills_meansHi <- expLogAndMeans(expressionGillsHi); expressionGills_meansLo <- expLogAndMeans(expressionGillsLo)
expressionLPJ_means <- expLogAndMeans(expressionLPJ); expressionLPJ_meansHi <- expLogAndMeans(expressionLPJHi); expressionLPJ_meansLo <- expLogAndMeans(expressionLPJLo)



pdf("~/CarolinGWAS/plots/expression_cacng.pdf",height=12,width=8)
par(mfrow=c(3,2))
barplot(expression_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[1],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(expressionLiver_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[2],ylim=c(0,7))
barplot(expressionOvary_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[3],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(expressionTestis_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[4],ylim=c(0,7))
barplot(expressionGills_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[5],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(expressionLPJ_means[cacng_gene_positions],names= cacng_genes$gene,las=2,main= tissues[6],ylim=c(0,7))
dev.off()


pdf("~/CarolinGWAS/plots/expression_cacng_behaviour_species_sameScale_hiLoComparison.pdf",height=12,width=8)
par(mfrow=c(3,2))
barplot(t(cbind(expression_meansHi[cacng_gene_positions],expression_meansLo[cacng_gene_positions])),beside=TRUE, names= cacng_genes$gene,las=2,main= tissues[1],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionLiver_meansHi[cacng_gene_positions],expressionLiver_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[2],ylim=c(0,7))
legend("topright",c("high exploration species","low exploration species"),col=c("grey10","grey70"),pch=15,cex=1.5)
barplot(t(cbind(expressionOvary_meansHi[cacng_gene_positions],expressionOvary_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[3],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionTestis_meansHi[cacng_gene_positions],expressionTestis_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[4],ylim=c(0,7))
barplot(t(cbind(expressionGills_meansHi[cacng_gene_positions],expressionGills_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[5],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionLPJ_meansHi[cacng_gene_positions],expressionLPJ_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[6],ylim=c(0,7))
dev.off()

pdf("~/CarolinGWAS/plots/expression_cacng_behaviour_species_sameScale_cacngSNPcomparison.pdf",height=12,width=8)
par(mfrow=c(3,2))
barplot(t(cbind(expression_meansHi[cacng_gene_positions],expression_meansLo[cacng_gene_positions])),beside=TRUE, names= cacng_genes$gene,las=2,main= tissues[1],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionLiver_meansHi[cacng_gene_positions],expressionLiver_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[2],ylim=c(0,7))
legend("topright",c("cacng5b SNP 0/0","cacng5b SNP 1/1"),col=c("grey10","grey70"),pch=15,cex=1.5)
barplot(t(cbind(expressionOvary_meansHi[cacng_gene_positions],expressionOvary_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[3],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionTestis_meansHi[cacng_gene_positions],expressionTestis_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[4],ylim=c(0,7))
barplot(t(cbind(expressionGills_meansHi[cacng_gene_positions],expressionGills_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[5],ylab="log2 expression (TMP)",ylim=c(0,7))
barplot(t(cbind(expressionLPJ_meansHi[cacng_gene_positions],expressionLPJ_meansLo[cacng_gene_positions])),beside=TRUE,names= cacng_genes$gene,las=2,main= tissues[6],ylim=c(0,7))
dev.off()





barplot(expression_meansHi[cacng_gene_positions])
