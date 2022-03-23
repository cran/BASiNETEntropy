#'@title Calculates the maximum entropy
#'@name maxentropy
#'
#'@description A function that calculates the maximum entropy
#'
#'@param histogram The histogram (used in 'training' function)
#'
#'@return Returns the maximum entropy
#'@author Murilo Montanini Breve
maxentropy <- function(histogram) {
	totalofpixels <- sum(histogram)
	maximum_entropy <- 0
	threshold <- NULL
	edgesname <- names(histogram)
	histogram <- unname(histogram)
	descendinghistogram <- sort(histogram, decreasing = TRUE)
	curveofentropy <- NULL
	for (t in seq_len(4095)) {
		P0 <- 0
		P1 <- 0
		for (i in seq_len(t)) {
			P0 <- P0 + descendinghistogram[i] / totalofpixels
		}
		for (i in (t + 1):4096) {
			P1 <- P1 + descendinghistogram[i] / totalofpixels
		}
		H0 <- 0
		H1 <- 0
		HT <- 0
		H0 <- H0 + entropy(P0)
		H1 <- H1 + entropy(P1)
		if (is.nan(H1) == TRUE)
			H1 <- 0
		if (is.nan(H0) == TRUE)
			H0 <- 0

		HT <- H0 + H1
		curveofentropy <- c(curveofentropy, HT)
		if (HT > maximum_entropy) {
			maximum_entropy <- HT
			threshold <- t
		}
		frequency <- descendinghistogram[threshold]
	}
	list <- list(maximum_entropy, threshold, frequency, curveofentropy)
	names(list) <-c("Max Entropy",
	                "Threshold",
	                "Edge frequency",
	                "Curve of Entropy")
	return(list)
}
