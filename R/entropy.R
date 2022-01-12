#'@title Calculates the entropy
#'@name entropy
#'
#'@description A function that calculates the entropy
#'
#'@param x The probabilities P0 and P1
#'
#'@return Returns the entropy
#'@author Murilo Montanini Breve
entropy <- function(x){
	entrop <- 0
	if (x > 0)
		entrop <- (-(x * log2(x)))
	if (is.nan(x) == TRUE)
		entrop <- 0
	return(entrop)
}
