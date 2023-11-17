



getGenesForGO <- function(GOdat,sigGenes,GOa) {
	genesWGO <- sigGenes[which(as.character(sigGenes) %in% genesInTerm(GOdat, GOa)[[1]])]
	genes <- getBM(c("entrezgene_id","description","external_gene_name","chromosome_name","start_position","end_position"),filters="entrezgene_id",values=genesWGO,mart = ensembl)
	#desc <- GenesLink[which(GenesLink[,7] %in% entrezAcc),]
	genes 
}


doOntology.CustomAnnot <- function(geneUniverse, sigGenes, GOmap, minNode=5) {
	testList <- factor(as.integer(unique(as.character(geneUniverse)) %in% unique(as.character(sigGenes))))
	names(testList) <- unique(geneUniverse)
	GOdata <- new("topGOdata", ontology = "MF", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	GOdataCC <- new("topGOdata", ontology = "CC", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	GOdataBP <- new("topGOdata", ontology = "BP", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	resultFisherWeigth <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	resultFisherCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
	allRes <- GenTable(GOdata, classic = resultFisher, weight = resultFisherWeigth, orderBy = "weight", ranksOf = "classic", topNodes = 20)
	resultFisherWeigthCC <- runTest(GOdataCC, algorithm = "weight", statistic = "fisher")
	allResCC <- GenTable(GOdataCC, classic = resultFisherCC, weight = resultFisherWeigthCC, orderBy = "weight", ranksOf = "classic", topNodes = 20)
	resultFisherBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
	resultFisherWeigthBP <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher")
	allResBP <- GenTable(GOdataBP, classic = resultFisherBP, weight = resultFisherWeigthBP, orderBy = "weight", ranksOf = "classic", topNodes = 40)
	list(allRes,allResCC,allResBP,GOdata,GOdataCC,GOdataBP)
}

doOntology.CustomAnnot.withResObjects <- function(geneUniverse, sigGenes, GOmap, minNode=5) {
	testList <- factor(as.integer(unique(as.character(geneUniverse)) %in% unique(as.character(sigGenes))))
	names(testList) <- unique(geneUniverse)
	GOdata <- new("topGOdata", ontology = "MF", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	GOdataCC <- new("topGOdata", ontology = "CC", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	GOdataBP <- new("topGOdata", ontology = "BP", allGenes = testList, annot = annFUN.gene2GO, gene2GO = GOmap, nodeSize= minNode)
	resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
	resultFisherWeigth <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
	resultFisherCC <- runTest(GOdataCC, algorithm = "classic", statistic = "fisher")
	allRes <- GenTable(GOdata, classic = resultFisher, weight = resultFisherWeigth, orderBy = "weight", ranksOf = "classic", topNodes = 20)
	resultFisherWeigthCC <- runTest(GOdataCC, algorithm = "weight", statistic = "fisher")
	allResCC <- GenTable(GOdataCC, classic = resultFisherCC, weight = resultFisherWeigthCC, orderBy = "weight", ranksOf = "classic", topNodes = 20)
	resultFisherBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
	resultFisherWeigthBP <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher")
	allResBP <- GenTable(GOdataBP, classic = resultFisherBP, weight = resultFisherWeigthBP, orderBy = "weight", ranksOf = "classic", topNodes = 40)
	list(allRes,allResCC,allResBP,GOdata,GOdataCC,GOdataBP, resultFisherWeigth, resultFisherWeigthCC, resultFisherWeigthBP)
}






findGenesNearSNPs <- function(SNPs,GenesGR,gap=50000) {
	SNPsStartEnd <- cbind(SNPs, SNPs$coord + 1); names(SNPsStartEnd) <- c(names(SNPs),"coordPlusOne")
	SNPsGR <- makeGRangesFromDataFrame(SNPsStartEnd, keep.extra.columns = TRUE,start.field="coord", end.field="coordPlusOne", seqnames.field="chr")
	ovlGenes <-  mergeByOverlaps(GenesGR, SNPsGR,maxgap=gap)
	#GO <- doOntology.CustomAnnot(geneUniverse, sigHDRs,AstCalGOmap,10)
	#GO[[7]] <- sigHDRs
	return(ovlGenes)
}



numToNCBI_onKocher <- function(num) {
	if (num == 1) { return("NC_031965"); }
	if (num == 2) { return("NC_031966"); }
	if (num == 3) { return("NC_031967"); }
	if (num == 4) { return("NC_031968"); }
	if (num == 5) { return("NC_031969"); }
	if (num == 6) { return("NC_031970"); }
	if (num == 7) { return("NC_031971"); }
	if (num == 8) { return("NC_031972"); }
	if (num == 9) { return("NC_031973"); }
	if (num == 10) { return("NC_031974"); }
	if (num == 11) { return("NC_031975"); }
	if (num == 12) { return("NC_031976"); }
	if (num == 13) { return("NC_031977"); }
	if (num == 14) { return("NC_031978"); }
	if (num == 15) { return("NC_031979"); }
	if (num == 16) { return("NC_031980"); }
	if (num == 17) { return("NC_031987"); }
	if (num == 18) { return("NC_031981"); }
	if (num == 19) { return("NC_031982"); }
	if (num == 20) { return("NC_031983"); }
	if (num == 21) { return("NC_031984"); }
	if (num == 22) { return("NC_031985"); }
	if (num == 23) { return("NC_031986"); }
	else { return(-1); }
}





LGtoNCBI_onKocher <- function(LG) {
	if (LG == "LG1") { return("NC_031965"); }
	if (LG == "LG2") { return("NC_031966"); }
	if (LG == "LG3a") { return("NC_031967"); }
	if (LG == "LG3b") { return("NC_031968"); }
	if (LG == "LG4") { return("NC_031969"); }
	if (LG == "LG5") { return("NC_031970"); }
	if (LG == "LG6") { return("NC_031971"); }
	if (LG == "LG7") { return("NC_031972"); }
	if (LG == "LG8") { return("NC_031973"); }
	if (LG == "LG9") { return("NC_031974"); }
	if (LG == "LG10") { return("NC_031975"); }
	if (LG == "LG11") { return("NC_031976"); }
	if (LG == "LG12") { return("NC_031977"); }
	if (LG == "LG13") { return("NC_031978"); }
	if (LG == "LG14") { return("NC_031979"); }
	if (LG == "LG15") { return("NC_031980"); }
	if (LG == "LG16") { return("NC_031987"); }
	if (LG == "LG17") { return("NC_031981"); }
	if (LG == "LG18") { return("NC_031982"); }
	if (LG == "LG19") { return("NC_031983"); }
	if (LG == "LG20") { return("NC_031984"); }
	if (LG == "LG22") { return("NC_031985"); }
	if (LG == "LG23") { return("NC_031986"); }
	else { return(LG); }
}


LGtoNCBIVec_onKocher <- function(LGvec) {
	NCBIvec <- numeric(0);
	for (i in 1:length(LGvec)) {
		NCBIvec <- c(NCBIvec,LGtoNCBI_onKocher(LGvec[i]))
	}
	return(NCBIvec)
}


# ------------ First declaring functions which we later use ---------------

lmp <- function (modelobject) {
    #if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

getGWASmiss <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    numNonMissing <- dim(modelobject$model)[1]
    if (numNonMissing == 45) {
    	return(TRUE)
    } else {
    	return(FALSE)
    }
}

doPGSL <- function(x) {
	ourDF <- data.frame(pheno= y$median_exploration,geno=as.numeric(x))
	rownames(ourDF) <- y$species_abb
	ourDF$species_abb = row.names(ourDF)
	compD = comparative.data(phy = tanTree2, data=ourDF, species_abb)
	if (length(unique(x)) == 1) { print(x); return(1.0); }
	else {
		fitmod= pgls( pheno ~ geno, compD)
		f <- summary(fitmod)$fstatistic
		p <- pf(f[1],f[2],f[3],lower.tail=F)
		attributes(p) <- NULL
		return(p)
	}
}


doPGSLexpression <- function(x,ourPheno,ourTree) {
	ourDF <- data.frame(pheno= ourPheno$median_exploration,geno=as.numeric(x))
	rownames(ourDF) <- ourPheno$species_abb
	ourDF$species_abb = row.names(ourDF)
	compD = comparative.data(phy = ourTree, data=ourDF, species_abb)
	if (length(unique(x)) == 1) { print(x); return(1.0); }
	else {
		fitmod= pgls( pheno ~ geno, compD)
		f <- summary(fitmod)$fstatistic
		p <- pf(f[1],f[2],f[3],lower.tail=F)
		attributes(p) <- NULL
		return(p)
	}
}




# Make a GO <-> gene list for the graph visualisation software
makeGOgeneCytoscape <- function() {
	xx <- as.list(org.Dr.egGO2EG); n <- names(xx)
	xy <- as.list(org.Dr.egGO2ALLEGS); n <- names(xy)
	xz <- as.list(GOTERM)
	allTerms <- numeric(0);
	allTermsSimple <- numeric(0);
	allTermsSimple3 <- numeric(0);
	for (i in 1:length(xz)) { 
		ei <- which(n == xz[[i]]@GOID);
		if (length(ei) >= 1) {
			entrez <-  paste(xy[[ei]],collapse="\t")
			d <- c(paste(toupper(xz[[i]]@Term),paste("GO",xz[[i]]@Ontology,sep=""),xz[[i]]@GOID,sep="%"),xz[[i]]@Term,entrez) 
			e <- c(xz[[i]]@GOID,xz[[i]]@Term,entrez)
			f <- c(xz[[i]]@Term,xz[[i]]@Definition,entrez)
			allTerms <- rbind(allTerms,d)
			allTermsSimple <- rbind(allTermsSimple,e)
			allTermsSimple3 <- rbind(allTermsSimple3,f);
		}
	}
	
	write.table(allTerms,"/Users/milanmalinsky/CarolinGWAS/EM_geneSetGWAS.gmt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(allTermsSimple,"/Users/milanmalinsky/CarolinGWAS/EM_geneSet_SimpleGWAS.gmt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	allTermsSimple2 <- cbind(allTermsSimple[,2],allTermsSimple[,2:3])
	write.table(allTermsSimple2,"/Users/milanmalinsky/CarolinGWAS/EM_geneSet_Simple2GWAS.gmt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(allTermsSimple3,"/Users/milanmalinsky/CarolinGWAS/EM_geneSet_Simple3GWAS.gmt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}



doOntologyOnPermutedSNPs <- function(numSNPs, allGenesGR, genesOnChroms, onGOmap) {
	snps <- permuteGWASsnps(numSNPs)
	near5kbGenes_snps <- findGenesNearSNPs(snps, allGenesGR, 5000)
	near5kbGenes_snps_IDs <- unique(near5kbGenes_snps$name2)
	GO_near5kbGenes_snps_5kb <- doOntology.CustomAnnot.withResObjects(genesOnChroms$name2, near5kbGenes_snps_IDs, onGOmap,10)
	return(GO_near5kbGenes_snps_5kb)
}




permuteGWASsnps <- function(numSNPs) {	
	permutedSNPs <- numeric(0)
	for (i in 1:numSNPs) {
		newChrNum <- sample(c(1:length(onKocherSizes[-1])),1, prob=onKocherSizes[-1]/sum(onKocherSizes))
		newChr <- numToNCBI_onKocher(newChrNum)
		maxLoc <- onKocherSizes[-1][newChrNum]
		newLoc <- sample(1:maxLoc,1)
		if (i == 1) {
			permutedSNPs <- c(newChr, newLoc)
		} else {
			permutedSNPs <- rbind(permutedSNPs,c(newChr, newLoc))
		}
	}
	permutedSNPsDF <- data.frame(chr=permutedSNPs[,1],coord=as.numeric(permutedSNPs[,2]))
	return(permutedSNPsDF)
}

permuteHDRs <- function(HDRs) {
	HDRlengths <- HDRs$nPos2-HDRs$nPos1
	sa <- t(sapply(HDRlengths, permuteHDR))
	newHDRs <- cbind(sa, HDRs$hdr, HDRs$val)
	colnames(newHDRs) <- colnames(HDRs)
	newHDRsDf <- as.data.frame(newHDRs)
	newHDRsDf$nPos1 <- as.numeric(as.character(newHDRsDf$nPos1))
	newHDRsDf$nPos2 <- as.numeric(as.character(newHDRsDf$nPos2))
	newHDRsDf$nChr...1. <- as.character(newHDRsDf$nChr...1.)
	return(newHDRsDf)
}

GOcorrectedForLengthGWAS <- function(allGenes,sigGenes) {
	AllGeneLengths <- allGenes$txEnd - allGenes$txStart
	DEgoSeq <- rep(0,length(allGenes$name2)); names(DEgoSeq) <- allGenes$name2
	DEgoSeq[which(names(DEgoSeq) %in% sigGenes)] <- 1;
	pwf <- nullp(DEgoSeq,bias.data= AllGeneLengths)
	go <- goseq(pwf,gene2cat=onGOmap,test.cats=c("GO:BP"))
	go <- go[which(go$ontology == "BP" & go$numInCat > 10),]
	goSig <- go[go$over_represented_pvalue<0.01,]
	goSigBHCorrected <- go[p.adjust(go$over_represented_pvalue, method="BH")<.05,]
	return(list(go[1:10,],goSig, goSigBHCorrected))
}


GOcorrectedForLengthGWAS_CC <- function(allGenes,sigGenes) {
	AllGeneLengths <- allGenes$txEnd - allGenes$txStart
	DEgoSeq <- rep(0,length(allGenes$name2)); names(DEgoSeq) <- allGenes$name2
	DEgoSeq[which(names(DEgoSeq) %in% sigGenes)] <- 1;
	pwf <- nullp(DEgoSeq,bias.data= AllGeneLengths)
	go <- goseq(pwf,gene2cat=onGOmap,test.cats=c("GO:CC"))
	go <- go[which(go$ontology == "CC" & go$numInCat > 10),]
	goSig <- go[go$over_represented_pvalue<0.01,]
	goSigBHCorrected <- go[p.adjust(go$over_represented_pvalue, method="BH")<.05,]
	return(list(go[1:10,],goSig, goSigBHCorrected))
}

GOcorrectedForLengthGWAS_MF <- function(allGenes,sigGenes) {
	AllGeneLengths <- allGenes$txEnd - allGenes$txStart
	DEgoSeq <- rep(0,length(allGenes$name2)); names(DEgoSeq) <- allGenes$name2
	DEgoSeq[which(names(DEgoSeq) %in% sigGenes)] <- 1;
	pwf <- nullp(DEgoSeq,bias.data= AllGeneLengths)
	go <- goseq(pwf,gene2cat=onGOmap,test.cats=c("GO:MF"))
	go <- go[which(go$ontology == "MF" & go$numInCat > 10),]
	goSig <- go[go$over_represented_pvalue<0.01,]
	goSigBHCorrected <- go[p.adjust(go$over_represented_pvalue, method="BH")<.05,]
	return(list(go[1:10,],goSig, goSigBHCorrected))
}


get5pUTRS <- function(x) {
	if (x[3] == "+") {
		return(c(x[2],x[4],x[6]))
	} else {
		return(c(x[2],x[7],x[5]))	
	}
}

get3pUTRS <- function(x) {
	if (x[3] == "+") {
		return(c(x[2],x[7],x[5]))
	} else {
		return(c(x[2],x[4],x[6]))	
	}
}

getUpstream5kbs <- function(x) {
	if (x[3] == "+") {
		return(c(x[2],as.numeric(x[4])-5000,as.numeric(x[4])-1))
	} else {
		return(c(x[2],as.numeric(x[5])+1,as.numeric(x[5])+5000))	
	}
}

getGeneBodies <- function(x) {
	if (x[3] == "+") {
		return(c(x[2],as.numeric(x[6]),as.numeric(x[7])+5000))
	} else {
		return(c(x[2],as.numeric(x[4])-5000,as.numeric(x[4])-1))	
	}
}

run_eWAS_GO_permutations <- function(ps_eWAS_LM,ofInterestIndex, onGOmap, sigGOIDsMF_eWAS, sigGOIDsCC_eWAS, sigGOIDsBP_eWAS) {
	for (i in 1:1000) {
		print(i)
		randomSampleOfGenes <- sample(names(ps_eWAS_LM),length(ofInterestIndex),replace=FALSE)
		GO_per_eWAS <- doOntology.CustomAnnot.withResObjects(names(ps_eWAS_LM), randomSampleOfGenes, onGOmap,10)
		allGoResMF_eWAS <- GenTable(GO_per_eWAS[[4]],weight=GO_per_eWAS[[7]], topNodes=length(usedGO(GO_per_eWAS[[4]])))
		MFpvals_eWAS <- allGoResMF_eWAS[which(allGoResMF_eWAS$GO.ID %in% sigGOIDsMF_eWAS),]
		MFpvals_eWAS <- MFpvals_eWAS[order(MFpvals_eWAS$GO.ID),]
		allGoResCC_eWAS <- GenTable(GO_per_eWAS[[5]],weight=GO_per_eWAS[[8]], topNodes=length(usedGO(GO_per_eWAS[[5]])))
		CCpvals_eWAS <- allGoResCC_eWAS[which(allGoResCC_eWAS$GO.ID %in% sigGOIDsCC_eWAS),]
		CCpvals_eWAS <- CCpvals_eWAS[order(CCpvals_eWAS$GO.ID),]
		allGoResBP_eWAS <- GenTable(GO_per_eWAS[[6]],weight=GO_per_eWAS[[9]], topNodes=length(usedGO(GO_per_eWAS[[6]])))
		BPpvals_eWAS <- allGoResBP_eWAS[which(allGoResBP_eWAS$GO.ID %in% sigGOIDsBP_eWAS),]
		BPpvals_eWAS <- BPpvals_eWAS[order(BPpvals_eWAS$GO.ID),]
	
		if (i == 1) {
			permutationResMF_eWAS <- cbind(sigGOIDsMF_eWAS[order(sigGOIDsMF_eWAS)], MFpvals_eWAS$Term, MFpvals_eWAS$weight)
			permutationResCC_eWAS <- cbind(sigGOIDsCC_eWAS[order(sigGOIDsCC_eWAS)], CCpvals_eWAS$Term, CCpvals_eWAS$weight)
			permutationResBP_eWAS <- cbind(sigGOIDsBP_eWAS[order(sigGOIDsBP_eWAS)], BPpvals_eWAS$Term, BPpvals_eWAS$weight)
		} else {
			permutationResMF_eWAS <- cbind(permutationResMF_eWAS, MFpvals_eWAS$weight)
			permutationResCC_eWAS <- cbind(permutationResCC_eWAS, CCpvals_eWAS$weight)
			permutationResBP_eWAS <- cbind(permutationResBP_eWAS, BPpvals_eWAS$weight)
		}
		
	} 
	
	permutationResMF_eWAS.onlyP <- matrix(as.numeric(permutationResMF_eWAS[,3:dim(permutationResMF_eWAS)[2]]),nrow=dim(permutationResMF_eWAS)[1])
	permutationResCC_eWAS.onlyP <- matrix(as.numeric(permutationResCC_eWAS[,3:dim(permutationResCC_eWAS)[2]]),nrow=dim(permutationResCC_eWAS)[1])
	permutationResBP_eWAS.onlyP <- matrix(as.numeric(permutationResBP_eWAS[,3:dim(permutationResBP_eWAS)[2]]),nrow=dim(permutationResBP_eWAS)[1])
	
	sigGOsMF_eWAS.ps <- as.numeric(sigGOsMF_eWAS[order(sigGOIDsMF_eWAS),]$weight)
	sigGOsCC_eWAS.ps <- as.numeric(sigGOsCC_eWAS[order(sigGOIDsCC_eWAS),]$weight)
	sigGOsBP_eWAS.ps <- as.numeric(sigGOsBP_eWAS[order(sigGOIDsBP_eWAS),]$weight)
	
	sigProportionsMF_eWAS <- numeric(0)
	sigProportionsCC_eWAS <- numeric(0)
	sigProportionsBP_eWAS <- numeric(0)
	for (i in 1:length(sigGOsMF_eWAS.ps)) {
		sigProportionsMF_eWAS <- c(sigProportionsMF_eWAS, length(which(permutationResMF_eWAS.onlyP[i,] <= sigGOsMF_eWAS.ps[i]))/length(permutationResMF_eWAS.onlyP[i,]))
		print(paste(permutationResMF_eWAS[i,2], sigProportionsMF_eWAS[i]))
	}
	
	for (i in 1:length(sigGOsCC_eWAS.ps)) {
		sigProportionsCC_eWAS <- c(sigProportionsCC_eWAS, length(which(permutationResCC_eWAS.onlyP[i,] <= sigGOsCC_eWAS.ps[i]))/length(permutationResCC_eWAS.onlyP[i,]))
		print(paste(permutationResCC_eWAS[i,2], sigProportionsCC_eWAS[i]))
	}
	
	for (i in 1:length(sigGOsBP_eWAS.ps)) {
		sigProportionsBP_eWAS <- c(sigProportionsBP_eWAS, length(which(permutationResBP_eWAS.onlyP[i,] <= sigGOsBP_eWAS.ps[i]))/length(permutationResBP_eWAS.onlyP[i,]))
		print(paste(permutationResBP_eWAS[i,2], sigProportionsBP_eWAS[i]))
	}
	
	write.table(permutationResMF_eWAS,"expression_GO_permutationResMF.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(permutationResCC_eWAS,"expression_GO_permutationResCC.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(permutationResBP_eWAS,"expression_GO_permutationResBP.txt",quote=F,sep="\t",row.names=F,col.names=F)
	
	write.table(cbind(permutationResMF_eWAS[,2], sigProportionsMF_eWAS),"expression_GO_permutationResMF_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cbind(permutationResCC_eWAS[,2], sigProportionsCC_eWAS),"expression_GO_permutationResCC_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cbind(permutationResBP_eWAS[,2], sigProportionsBP_eWAS),"expression_GO_permutationResBP_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
}



### ------------------------------------- Section on permutations ----------------------------------------------------
# --- This needed to be done only once - the results are witten into text files at the end of this section 
doGOpermutations <- function() {
	for (i in 1:1000) {
		print(i)
		GO_per <- doOntologyOnPermutedSNPs(dim(topLMsnpsPermissivePGLS)[1],allGenesGR,genesOnChroms, onGOmap)
		allGoResMF <- GenTable(GO_per[[4]],weight=GO_per[[7]], topNodes=length(usedGO(GO_per[[4]])))
		MFpvals <- allGoResMF[which(allGoResMF$GO.ID %in% sigGOIDsMF),]
		MFpvals <- MFpvals[order(MFpvals$GO.ID),]
		allGoResCC <- GenTable(GO_per[[5]],weight=GO_per[[8]], topNodes=length(usedGO(GO_per[[5]])))
		CCpvals <- allGoResCC[which(allGoResCC$GO.ID %in% sigGOIDsCC),]
		CCpvals <- CCpvals[order(CCpvals$GO.ID),]
		allGoResBP <- GenTable(GO_per[[6]],weight=GO_per[[9]], topNodes=length(usedGO(GO_per[[6]])))
		BPpvals <- allGoResBP[which(allGoResBP$GO.ID %in% sigGOIDsBP),]
		BPpvals <- BPpvals[order(BPpvals$GO.ID),]
		
		if (i == 1) {
			permutationResMF <- cbind(sigGOIDsMF[order(sigGOIDsMF)], MFpvals$weight, MFpvals$Term)
			permutationResCC <- cbind(sigGOIDsCC[order(sigGOIDsCC)], CCpvals$weight, CCpvals$Term)
			permutationResBP <- cbind(sigGOIDsBP[order(sigGOIDsBP)], BPpvals$weight, BPpvals$Term)
		} else {
			permutationResMF <- cbind(permutationResMF, MFpvals$weight)
			permutationResCC <- cbind(permutationResCC, CCpvals$weight)
			permutationResBP <- cbind(permutationResBP, BPpvals$weight)
		}
		
	} 
	permutationResMF.onlyP <- matrix(as.numeric(permutationResMF[,c(2,4:dim(permutationResMF)[2])]),nrow=dim(permutationResMF)[1])
	permutationResCC.onlyP <- matrix(as.numeric(permutationResCC[,c(2,4:dim(permutationResCC)[2])]),nrow=dim(permutationResCC)[1])
	permutationResBP.onlyP <- matrix(as.numeric(permutationResBP[,c(2,4:dim(permutationResBP)[2])]),nrow=dim(permutationResBP)[1])
	
	sigGOsMF.ps <- as.numeric(sigGOsMF[order(sigGOIDsMF),]$weight)
	sigGOsCC.ps <- as.numeric(sigGOsCC[order(sigGOIDsCC),]$weight)
	sigGOsBP.ps <- as.numeric(sigGOsBP[order(sigGOIDsBP),]$weight)
	
	
	sigProportionsMF <- numeric(0)
	sigProportionsCC <- numeric(0)
	sigProportionsBP <- numeric(0)
	for (i in 1:length(sigGOsMF.ps)) {
		sigProportionsMF <- c(sigProportionsMF, length(which(permutationResMF.onlyP[i,] <= sigGOsMF.ps[i]))/length(permutationResMF.onlyP[i,]))
		print(paste(permutationResMF[i,3], sigProportionsMF[i]))
	}
	
	for (i in 1:length(sigGOsCC.ps)) {
		sigProportionsCC <- c(sigProportionsCC, length(which(permutationResCC.onlyP[i,] <= sigGOsCC.ps[i]))/length(permutationResCC.onlyP[i,]))
		print(paste(permutationResCC[i,3], sigProportionsCC[i]))
	}
	
	for (i in 1:length(sigGOsBP.ps)) {
		sigProportionsBP <- c(sigProportionsBP, length(which(permutationResBP.onlyP[i,] <= sigGOsBP.ps[i]))/length(permutationResBP.onlyP[i,]))
		print(paste(permutationResBP[i,3], sigProportionsBP[i]))
	}
	
	write.table(permutationResMF[,c(1,3,2,4:1002)],"permutationResMF.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(permutationResCC[,c(1,3,2,4:1002)],"permutationResCC.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(permutationResBP[,c(1,3,2,4:1002)],"permutationResBP.txt",quote=F,sep="\t",row.names=F,col.names=F)
	
	write.table(cbind(permutationResMF[,3], sigProportionsMF),"permutationResMF_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cbind(permutationResCC[,3], sigProportionsCC),"permutationResCC_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
	write.table(cbind(permutationResBP[,3], sigProportionsBP),"permutationResBP_means.txt",quote=F,sep="\t",row.names=F,col.names=F)
}


