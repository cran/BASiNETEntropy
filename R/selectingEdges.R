#'@title Selects the edges of the adjacency matrix
#'@name selectingEdges
#'
#'@description A function that selects the edges of the adjacency matrix
#'
#'@param MAX The maximum entropy
#'@param data The adjacency matrix
#'
#'@return Returns the selected edges of the adjacency matrix
#'@author Murilo Montanini Breve
selectingEdges <- function(MAX, data) {
	for (u in seq_along(data)) {
		if (data[u] > MAX[[3]]) {
			data[u] <- 1
		} else
			data[u] <- 0
	}
	return(data)
}
