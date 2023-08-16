#'@title Performs the classification methodology using complex network and
#'entropy theories
#'@name classify
#'
#'@description Given three or two distinct data sets, one of mRNA, one
#'of lncRNA and one of sncRNA.
#'The classification of the data is done from the structure of the networks
#' formed by the sequences, that is filtered by an entropy methodology.
#'After this is done, the classification starts.
#'
#'@param mRNA Directory where the file .FASTA lies with the mRNA sequences
#'@param lncRNA Directory where the file .FASTA lies with the lncRNA sequences
#'@param sncRNA Directory where the file .FASTA lies with the sncRNA sequences
#' (optional)
#'@param trainingResult The result of the training, (three or two matrices)
#'@param save_dataframe save when set, this parameter saves a .csv file with
#'the features in the current directory. No file is created by default.
#'@param save_model save when set, this parameter saves a .rds file with
#'the model in the current directory. No file is created by default.
#'@param predict_with_model predict the input sequences with the previously
#'generated model.
#'
#'@return Results
#'
#'@author Murilo Montanini Breve
#'
#'@examples
#'library(BASiNETEntropy)
#'arqSeqMRNA <- system.file("extdata", "mRNA.fasta",package = "BASiNETEntropy")
#'arqSeqLNCRNA <- system.file("extdata", "ncRNA.fasta", package = "BASiNETEntropy")
#'load(system.file("extdata", "trainingResult.RData", package = "BASiNETEntropy"))
#'r_classify <- classify(mRNA=arqSeqMRNA, lncRNA=arqSeqLNCRNA, trainingResult = trainingResult)

#'@importFrom Biostrings readBStringSet
#'@import igraph
#'@import randomForest
#'@importFrom graphics abline axis legend text
#'@importFrom stats sd predict
#'@importFrom utils write.csv2
#'@export
classify <- function(mRNA,
                     lncRNA,
                     sncRNA = NULL,
                     trainingResult,
                     save_dataframe = NULL,
                     save_model = NULL,
                     predict_with_model = NULL) {
  if (missing(trainingResult))
    trainingResult <- training(mRNA, lncRNA, sncRNA)

  vectorm <- NULL
  vectorlnc <- NULL
  vectorsnc <- NULL
  edgeslistmrna <- NULL
  edgeslistlncrna <- NULL
  edgeslistsncrna <- NULL
  MRNA <- readBStringSet(mRNA)
  LNCRNA <- readBStringSet(lncRNA)
  if (length(sncRNA)) {
    SNCRNA <- readBStringSet(sncRNA)
  }
  for (t in 1:3) {
    if (t == 1) {
      message("[INFO] Classifying mRNA:")
      seq <- c(MRNA)
    }
    if (t == 2) {
      message("[INFO] Classifying lncRNA:")
      seq <- c(LNCRNA)
    }
    if (t == 3) {
      if (length(sncRNA)) {
        message("[INFO] Classifying sncRNA:")
        seq <- c(SNCRNA)
      } else
        break
    }
    for (u in seq_along(seq)) {
      sequence <- strsplit(toString(seq[u]), split = '')
      sequence <- sequence[[1]]
      if (t == 1)
        edgeslistmrna[[u]] <- createedges(sequence)
      if (t == 2)
        edgeslistlncrna[[u]] <- createedges(sequence)
      if (t == 3)
        edgeslistsncrna[[u]] <- createedges(sequence)
      message(u)
    }
  }

  message("[INFO] Filtering the graphs")

  sequenciam <- filtering(trainingResult[[1]], edgeslistmrna)
  sequencial <- filtering(trainingResult[[2]], edgeslistlncrna)
  if (length(sncRNA)) {
    sequencias <- filtering(trainingResult[[3]], edgeslistsncrna)
  }

  message("[INFO] Extracting measurements from graphs")
  measures <- NULL

  for(q in seq_along(sequenciam)){
    net <- graph_from_adjacency_matrix(sequenciam[[q]], mode = c("undirected"))
    measures<-c(measures,average.path.length(net,directed = FALSE,
                                             unconnected=FALSE))
    measures<-c(measures,transitivity(net,
                                      type=c("undirected"),
                                      vids=NULL,
                                      weights=NULL,
                                      isolates=c("NaN","zero")))
    measures<-c(measures,mean(degree(net, v=V(net),
                                     normalized=FALSE)))
    measures<-c(measures,assortativity_degree(net,
                                              directed = FALSE))
    measures<-c(measures,mean(betweenness(net,
                                          v = V(net),
                                          directed = FALSE,
                                          weights = NULL,
                                          normalized = FALSE)))
    measures<-c(measures,sd(degree(net,
                                   v=V(net),
                                   normalized = FALSE),
                            na.rm = FALSE))
    measures<-c(measures,which.max(degree(net,
                                          v=V(net),
                                          normalized=FALSE)))
    measures<-c(measures,which.min(degree(net,
                                          v = V(net),
                                          normalized=FALSE)))
    measures<-c(measures,(count_motifs(net, size = 3)))
    measures<-c(measures,(count_motifs(net, size = 4)))
  }

  for(q in seq_along(sequencial)){
    net<-graph_from_adjacency_matrix(sequencial[[q]],
                                     mode = c("undirected"))
    measures<-c(measures,average.path.length(net,
                                             directed = FALSE,
                                             unconnected=FALSE))
    measures<-c(measures,transitivity(net,
                                      type = c("undirected"),
                                      vids=NULL,
                                      weights=NULL,
                                      isolates=c("NaN","zero")))
    measures<-c(measures,mean(degree(net,
                                     v = V(net),
                                     normalized=FALSE)))
    measures<-c(measures,assortativity_degree(net,
                                              directed = FALSE))
    measures<-c(measures,mean(betweenness(net,
                                          v = V(net),
                                          directed = FALSE,
                                          weights = NULL,
                                          normalized = FALSE)))
    measures<-c(measures,sd(degree(net,
                                   v = V(net),
                                   normalized= FALSE),
                            na.rm = FALSE))
    measures<-c(measures,which.max(degree(net,
                                          v=V(net),
                                          normalized=FALSE)))
    measures<-c(measures,which.min(degree(net,
                                          v = V(net),
                                          normalized=FALSE)))
    measures<-c(measures,(count_motifs(net, size = 3)))
    measures<-c(measures,(count_motifs(net, size = 4)))
  }

  if(length(sncRNA)) {
    for(q in seq_along(sequencias)){
      net<-graph_from_adjacency_matrix(sequencias[[q]],
                                       mode = c("undirected"))
      measures<-c(measures,average.path.length(net,
                                               directed=FALSE,
                                               unconnected=FALSE))
      measures<-c(measures,transitivity(net,
                                        type = c("undirected"),
                                        vids=NULL,
                                        weights=NULL,
                                        isolates=c("NaN","zero")))
      measures<-c(measures,mean(degree(net,
                                       v = V(net),
                                       normalized=FALSE)))
      measures<-c(measures,assortativity_degree(net,
                                                directed = FALSE))
      measures<-c(measures,mean(betweenness(net,
                                            v = V(net),
                                            directed = FALSE,
                                            weights = NULL,
                                            normalized = FALSE)))
      measures<-c(measures,sd(degree(net,
                                     v = V(net),
                                     normalized= FALSE),
                              na.rm = FALSE))
      measures<-c(measures,which.max(degree(net,
                                            v = V(net),
                                            normalized=FALSE)))
      measures<-c(measures,which.min(degree(net,
                                            v = V(net),
                                            normalized=FALSE)))
      measures<-c(measures,(count_motifs(net, size = 3)))
      measures<-c(measures,(count_motifs(net, size = 4)))
    }}

  message("[INFO] Building the dataframes")

  if (length(sncRNA)) {
    data <- creatingDataframe(measures,
                              tamM = length(sequenciam),
                              tamLNC = length(sequencial),
                              tamSNC = length(sequencias))
  } else
    data <- creatingDataframe(measures,
                              tamM = length(sequenciam),
                              tamLNC = length(sequencial))
  DF <- data[,-11]
  DF[is.na(DF)] <- 0

  if (!missing(save_dataframe)) {
    write.csv2(data, file = "feature_matrix.csv")
    message("feature_matrix.csv file generated in the current R directory")
  }

  if (missing(predict_with_model)) {
    message("[INFO] Sorting with Randomforest")
    rf <- randomForest(DF, as.factor(data[, 11]))
    print(rf)

    if (!missing(save_model)) {
      save(rf, file = paste(save_model, ".Rda", sep = ""))
      message(
        paste(predict_with_model, ".Rda", sep = ""),
        " file generated in the current R directory"
      )
    }

    return(rf)

  } else{

    message("[INFO] Predicting with the input model")
    load(file = paste(predict_with_model, ".Rda", sep = ""))
    rf.pred <-  predict(rf, DF)
    print(rf.pred)

    return(rf.pred)

  }

}
