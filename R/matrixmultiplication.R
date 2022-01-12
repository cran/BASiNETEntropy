#'@title Compares the matrices
#'@name matrixmultiplication
#'
#'@description A function that compares the matrices 'trainingResult' and
#'the adjacency matrix to produce a filtered adjacency matrix.
#'
#'@param data Adjacency matrix
#'@param histodata 'trainingResult' data
#'
#'@return Returns the filtered adjacency matrix
#'@author Murilo Montanini Breve
matrixmultiplication <- function(data, histodata) {
	ordereddata <- data[order(rownames(data)),
	                    order(colnames(data))]
	orderedhistogram <-	histodata[order(rownames(histodata)),
	                              order(colnames(histodata))]
	datanames <- rownames(ordereddata)
	histonames <- rownames(orderedhistogram)
	selectednames <- intersect(histonames, datanames)
	ordereddata <- ordereddata[selectednames, selectednames]
	selectedmatrix <-	matrix(nrow = length(ordereddata[, 1]),
	                         ncol = length(ordereddata[, 1]))

	orderedhistogram <- orderedhistogram[selectednames, selectednames]
	for (i in seq_along(ordereddata[, 1])) {
		for (j in seq_along(ordereddata[, 1])) {
			selectedmatrix[i, j] <- ordereddata[i, j] * orderedhistogram[i, j]
		}
	}
	rownames(selectedmatrix) <- rownames(ordereddata)
	colnames(selectedmatrix) <- rownames(ordereddata)
	return(selectedmatrix)
}
