#'@title Creates an untargeted graph from a biological sequence
#'@name createedges
#'
#'@description A function that from a biological sequence generates a graph
#' not addressed having as words vertices, this being able to have its size
#' parameter set by the' word 'parameter. The connections between words depend
#' of the' step 'parameter that indicates the next connection to be formed
#'
#'@param step It is the integer parameter that decides the step that will
#'be taken to make a new connection
#'@param word This integer parameter decides the size of the word that will
#'be formed
#'@param sequence It is a vector that represents the sequence
#'
#'@return Returns the array used to creates the edge list
#'@author Murilo Montanini Breve
createedges <- function(sequence, word = 3, step = 1) {
	aux <- ""
	index <- 1
	position <- 0
	cont <- length(sequence)
	comma <- 0
	x <- 0
	k <- 1
	vector <- c()
	while ((index - 1 + (word * 2)) < cont){
		while (x < word){
			aux <- paste(aux, sequence[index], sep = "")
			x <- x + 1
			index <- index + 1
		}
		vector <- c(vector, aux)
		aux <- ""
		x <- 0
		while (x < word){
			aux <- paste(aux, sequence[index], sep = "")
			x <- x + 1
			index <- index + 1
		}
		vector <- c(vector, aux)
		aux <- ""
		x <- 0
		position <- position + step
		index <- position + 1
	}
	return(vector)
}
