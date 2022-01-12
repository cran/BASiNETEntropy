#'@title Filters the edges
#'@name filtering
#'
#'@description A function that filters the edges after the maximum entropy is
#' obtained
#'
#'@param edgestoselect The selected edges
#'@param edgestofilter The edges used to filter
#'
#'@return Returns the filtered edges
#'@author Murilo Montanini Breve
#'@import igraph
filtering <- function(edgestoselect, edgestofilter) {
	filtered <- NULL
	for (t in seq_along(edgestofilter)) {
		net <- graph(edges = edgestofilter[[t]], directed = FALSE)
		matrix <- as_adjacency_matrix(net)
		data <- as.matrix(matrix)
		filtered[[t]] <- matrixmultiplication(data, edgestoselect)
	}
	return(filtered)
}
