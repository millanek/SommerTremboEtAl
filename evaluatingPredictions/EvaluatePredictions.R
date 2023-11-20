########################################################################
### Script by Milan Malinsky
### Last edit: 20th Nov 2023
### Exploring the gene expression of cacng genes across tissues

################################################################################
##  input data        				      ######################################
################################################################################

setwd("~/CarolinGWAS/SommerTremboEtAl/evaluatingPredictions/") 

y <- read.table("../TestingForAssociation/exploratoryBehaviorMedians.txt",header=T)
predictions <- c(0.105, 0.162, 0.195, 0.125, 0.101, 0.095, 0.366, 0.378)
observations <- c(0.1402715, 0.03280543, 0.02714932, 0.09615385, 0.04638009, 0.04524887, 0.4864253, 0.6945701)
predictionsBimodal <- c(0, 0, 0, 0, 0, 0, 1, 1)

################################################################################
## 		Accuracy of the continuous behavioral predictions   ####################
################################################################################

deltaData <- mean(abs(predictions - observations))  # The difference between predictions and real data

# A million random subsamples of eight scores from among the 56 initial species
b <- y$median_exploration
x <- t(replicate(1000000,sample(b,8)))
diffsReplicates <- abs(sweep(x,2, predictions))
deltasReplicates <- apply(diffsReplicates, 1, mean)

# Fig. S16
h <- hist(deltasReplicates,breaks=100,xlab="Average prediction error",freq=T,yaxt='n',ylab="Proportion of random samples",main="")
abline(v= deltaData,col="red",lwd=2)
axis(2,at=c(0,10000,20000,30000,40000,50000),labels=c("0","0.01","0.02","0.03","0.04","0.05"))
legend("topright","Prediction error for real\nempirical values = 0.113",lty=1,col="red",lwd=2,bty = "n")
edcfReplicates <- ecdf(deltasReplicates); edcfReplicates(deltaData)

################################################################################
## 		Accuracy of binary predictions  					    ####################
################################################################################

# Threshold - bimodal predictions 
threshold <- 0.3535068          # mean between the highest C/C (Limdar) and the lowest T/T species (Intloo) at the top SNP
b[ b<threshold ] <- 0; b[ b>=threshold ] <- 1; # Make the behavioral values bimodal

# 10000 random subsamples of eight bimodal values from among the 56 initial species
x <- t(replicate(10000,sample(b,8)))
diffsReplicates <- abs(sweep(x,2, predictionsBimodal))
deltasReplicates <- apply(diffsReplicates, 1, sum)

# Fig. S15
sumsDistances <- table(deltasReplicates)
sumsDistances <- sumsDistances / length(deltasReplicates)
names(sumsDistances) <- c("8","7","6","5","4","3","2","1")
barplot(rev(sumsDistances),xlab = "Number of correct precitions",ylab="Probability at random")  
length(which(deltasReplicates == 0))/length(deltasReplicates)


