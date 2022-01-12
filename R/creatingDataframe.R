#'@title Creates a feature matrix using complex network topological measures
#'@name creatingDataframe
#'
#'@description A function that from the complex network topological measures
#'create the feature matrix.
#'
#'@param measures The complex network topological measures
#'@param tamM mRNA sequence size
#'@param tamLNC lncRNA sequence size
#'@param tamSNC snRNA sequence size
#'
#'@return Returns the feature matrix in scale 0-1
#'@author Murilo Montanini Breve
creatingDataframe <- function(measures, tamM, tamLNC, tamSNC) {
    if (missing(tamSNC)) {
        Specie <- NULL
        for (i in 1:tamM)
            Specie <- c(Specie, "mRNA")
        for (i in (tamM + 1):(tamM + tamLNC))
            Specie <- c(Specie, "lncRNA")
        dataframe <- NULL
        dataframe <- matrix(data = measures,
                            nrow = (tamM + tamLNC),
                            ncol = 10,
                            byrow = TRUE | TRUE | FALSE | FALSE)

        dataframe <- preprocessing(dataframe, tamM, tamLNC)
        Species <- matrix(data = Specie,
                          nrow = (tamM + tamLNC),
                          ncol = 1,
                          byrow = TRUE | TRUE | FALSE | FALSE)
        dataframe <- cbind(dataframe, Species)
        rownames(dataframe) <- Specie
        colnames(dataframe) <- c("ASPL",
                                 "CC",
                                 "DEG",
                                 "ASS",
                                 "BET",
                                 "SD",
                                 "MAX",
                                 "MIN",
                                 "MT3",
                                 "MT4",
                                 "CLASS")
        DF = as.data.frame(dataframe)
    }
    else{
        Specie <- NULL
        for (i in 1:tamM)
            Specie <- c(Specie, "mRNA")
        for (i in (tamM + 1):(tamM + tamLNC))
            Specie <- c(Specie, "lncRNA")
        for (i in (tamM + tamLNC + 1):(tamM + tamLNC + tamSNC))
            Specie <- c(Specie, "sncRNA")
        dataframe <- NULL
        dataframe <- matrix(data = measures,
                            nrow = (tamM + tamLNC + tamSNC),
                            ncol = 10,
                            byrow = TRUE | TRUE | FALSE | FALSE)
        dataframe <- preprocessing(dataframe, tamM, tamLNC, tamSNC)
        Species <- matrix(data = Specie,
                          nrow = (tamM + tamLNC + tamSNC),
                          ncol = 1,
                          byrow = TRUE | TRUE | FALSE | FALSE)
        dataframe <- cbind(dataframe, Species)
        rownames(dataframe) <- Specie
        colnames(dataframe) <- c("ASPL",
                                 "CC",
                                 "DEG",
                                 "ASS",
                                 "BET",
                                 "SD",
                                 "MAX",
                                 "MIN",
                                 "MT3",
                                 "MT4",
                                 "CLASS")
        DF = as.data.frame(dataframe)
    }
    return(DF)
}
