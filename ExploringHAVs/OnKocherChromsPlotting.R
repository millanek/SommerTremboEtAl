

posToLinearPos <- function(chr,pos) {
	absPos <- pos + onKocherSizesSums[chr] + (chr-1)* 5000000
	return(absPos)
}
chrMids <- function(onKocherSizes, onKocherSizesSums,gap) {
	mids <- numeric(0);
	for (i in 1:23) {
		mids[i] <- onKocherSizesSums[i] + (onKocherSizes[i+1]/2) + (i-1)* 5000000
	}
	return(mids)
}

onKocherSizes <- c(0,38372991,35256741,14041792,54508961,38038224,34628617,44571662,62059223, 30802437, 27519051, 32426571, 36466354, 41232431, 32337344, 39264731, 36154882, 43860769, 40919683, 37007722, 31245232, 36767035, 37011614,44097196)
onKocherSizesSums <- cumsum(onKocherSizes)
onKocherMids <- chrMids(onKocherSizes, onKocherSizesSums,5000000)