#'@title Trains the algorithm to select the edges that maximize the entropy
#'@name training
#'
#'@description A function that trains the algorithm to select the edges that
#'maximize the entropy
#'
#'@param mRNA Directory where the file .FASTA lies with the mRNA sequences
#'@param lncRNA Directory where the file .FASTA lies with the lncRNA sequences
#'@param sncRNA Directory where the file .FASTA lies with the sncRNA sequences
#' (optional)
#'
#'@return Returns the edge lists and the 'curveofentropy' function inputs
#'@author Murilo Montanini Breve
#'
#'@importFrom Biostrings readBStringSet
#'@import igraph
#'@import randomForest
#'@importFrom graphics abline axis legend text
#'@importFrom stats sd
#'@importFrom utils write.csv2
training <- function(mRNA, lncRNA, sncRNA = NULL) {
	MRNA <- readBStringSet(mRNA)
	LNCRNA <- readBStringSet(lncRNA)
	if (length(sncRNA)) {
		SNCRNA <- readBStringSet(sncRNA)
	}
	vectorm <- NULL
	vectorlnc <- NULL
	vectorsnc <- NULL
	edgeslistmrna <- NULL
	edgeslistlncrna <- NULL
	edgeslistsncrna <- NULL
	mRNAmatrixEdges <- NULL
	lncRNAmatrixEdges <- NULL
	sncRNAmatrixEdges <- NULL
	MAXsnc <- NULL

	for (t in seq_len(3)) {
		if (t == 1) {
			message("[INFO] Analyzing mRNA:")
			seq <- c(MRNA)
		}
		if (t == 2) {
			message("[INFO] Analyzing lncRNA:")
			seq <- c(LNCRNA)
		}
		if (t == 3) {
			if (length(sncRNA)) {
				message("[INFO] Analyzing sncRNA:")
				seq <- c(SNCRNA)
			} else
				break
		}
		for (u in seq_along(seq)) {
			sequence <- strsplit(toString(seq[u]), split = '')
			sequence <- sequence[[1]]
			aux <- ""
			index <- 1
			position <- 0
			cont <- length(sequence)
			comma <- 0
			x <- 0
			k <- 1
			vector <- c()
			if (t == 1)
				edgeslistmrna[[u]] <- createedges(sequence)
			if (t == 2)
				edgeslistlncrna[[u]] <- createedges(sequence)
			if (t == 3)
				edgeslistsncrna[[u]] <- createedges(sequence)
			message(u)
		}
	}
	for (f in seq_along(edgeslistmrna)) {
		vectorm <- c(vectorm, edgeslistmrna[[f]])
	}
	for (f in seq_along(edgeslistlncrna)) {
		vectorlnc <- c(vectorlnc, edgeslistlncrna[[f]])
	}
	if (length(sncRNA)) {
		for (f in seq_along(edgeslistsncrna)) {
			vectorsnc <- c(vectorsnc, edgeslistsncrna[[f]])
		}
	}
	netm <- graph(edges = vectorm, directed = FALSE)
	netl <- graph(edges = vectorlnc, directed = FALSE)
	matrizm <- as_adjacency_matrix(netm)
	datam <- as.matrix(matrizm)
	matrizl <- as_adjacency_matrix(netl)
	datal <- as.matrix(matrizl)
	datam <- datam[order(rownames(datam)), order(colnames(datam))]
	datal <- datal[order(rownames(datal)), order(colnames(datal))]
	message("[INFO] Analyzing entropy")
	MAXm <- maxentropy(datam)
	MAXlnc <- maxentropy(datal)

	if (length(sncRNA)) {
		nets <- graph(edges = vectorsnc, directed = FALSE)
		matrizs <- as_adjacency_matrix(nets)
		datas <- as.matrix(matrizs)
		MAXsnc <- maxentropy(datas)
		sncRNAmatrixEdges <- selectingEdges(MAXsnc, datas)
	}

	message("[INFO] Selecting the edges by the maximum entropy method")
	mRNAmatrixEdges <- selectingEdges(MAXm, datam)
	lncRNAmatrixEdges <- selectingEdges(MAXlnc, datal)
	listMatrix <- list(mRNAmatrixEdges,
	                   lncRNAmatrixEdges,
	                   sncRNAmatrixEdges,
	                   MAXm,
	                   MAXlnc,
	                   MAXsnc)
	return(listMatrix)
}
