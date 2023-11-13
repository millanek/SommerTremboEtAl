########################################################################
### Script by Milan Malinsky and Carolin Sommer-Trembo
### Last edit: 17th June 2021
### This is an executable R script for running GWAS with correction for phylogeny (but not gene-flow) using the PGLS method
### GWAS means associating species-level allele frequencies with the behaviour median values for those species
### Takes some time (definitely many hours) to run for each chromosome
### Example to execute this script:
# Rscript PGLSrun.R NC_031969_AlleleFrequencies_fromProbabilities_n100000.txt.gz 

########################################################################
#######   Load packages  ######################################
########################################################################

library(caper)

########################################################################
#######   Defining functions that are used below  ######################
########################################################################

# Obviously - runs the PGLS
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

########################################################################
##  The script to be executed                 ##########################
########################################################################

############
### Read and process the Allele Frequencies from the file supplied by the command line argument
############

args = commandArgs(trailingOnly=TRUE)
AFfile <- args[1]
AF <- read.table(AFfile,header=T)
# -1 stands originally for missing data; then subset only the Allele Frequency columns
AF[AF==-1] <- NA; AFfull <- AF; AF <- AF[,5:60]
# Scale the actual Allele Frequncies to be on the interval [-1,1]  
AFscaled <- 2*((AF - min(AF,na.rm=T))/(max(AF,na.rm=T)-min(AF,na.rm=T))) - 1

############
### Read phenotype vector  
############ 
y <- read.table("exploratoryBehaviorMedians.txt",header=T)
y <- y[match(colnames(AFscaled), y$species_abb),]

##########
### Remove SNPs with missing data (i.e. those where the allele frequency could not be calculated)
#########
N <- apply(AFscaled,1, function(x) length(which(is.na(x == T))))        # This calculates the number of missing data per SNP
AFscaled.noMiss <- AFscaled[which(N == 0),];                            # This subsets the AF data to remove SNPs with missing data 

##########
# The way things are done here, we actually still need to run the standard LM GWAS: 
# y = X_a*beta_a + epsilon 
# in order to find some SNPs for which because of missing data we can't get the slope of the lm
# and then we exclude these SNPs from the PGLS runs below
##########
GWASvector.noMiss <- apply(AFscaled.noMiss, 1, function(x) lm(y$median_exploration ~ x))
naVector <- sapply(GWASvector.noMiss, function(x) is.na(x$coefficients[2]));    ### this object is used below in the PGLS to exclude the SNPs 

##########
### Here we do the actual PGLS!!!!!!! takes many hours
### We load the phenotype vector again and process it slightly differently for this purpose - this is not optimal, but it works so let's not mess with it for no good reason
#########
tanTree <- read.tree("TanganyikaSpeciesTree_b1.tre"); 
exploMeasurement <- read.table("exploratoryBehaviorMedians.txt",header=T); exploMeasurement <- exploMeasurement[-1,]
tanTree2 <- drop.tip(tanTree, tanTree$tip.label[! tanTree$tip.label %in% exploMeasurement$species_abb]);
pglsGWASvector <- apply(AFscaled.noMiss[-which(naVector == TRUE),], 1, doPGSL)
pglsGWASvector.withLoc <- cbind(AFfull[names(pglsGWASvector),1:2],pglsGWASvector)
write.table(pglsGWASvector.withLoc,file=paste0(substr(AFfile, 1, nchar(AFfile)-7),"_PGLSpiVals_corrected.txt"),quote=F,sep="\t")
