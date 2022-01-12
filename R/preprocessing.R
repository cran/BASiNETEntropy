#'@title Rescales the results between values from 0 to 1
#'@name preprocessing
#'
#'@description Given the results the data is rescaled for values between
#'0 and 1, so that the length of the sequences does not influence the results.
#'The rescaling of the sequences are made separately
#'
#'@param datah Array with results numerics
#'@param tamM Integer number of mRNA sequences
#'@param tamLNC Integer number of lncRNA sequences
#'@param tamSNC Integer number of sncRNA sequences
#'
#'@return Returns the array with the rescaled values
#'@author Murilo Montanini Breve
preprocessing <- function(datah, tamM, tamLNC, tamSNC) {
	range1 <- c()
	range2 <- c()
	range3 <- c()
	mini1 <- c()
	maxi1 <- c()
	mini2 <- c()
	mini3 <- c()
	maxi2 <- c()
	maxi3 <- c()

	if (missing(tamSNC)) {
		for (u in seq_len(10)) {
			range1 <- c(range1, range(datah[1:tamM, u]))
		}
		for (u in seq_len(10)) {
			range2 <- c(range2, range(datah[(tamM + 1):(tamM + tamLNC), u]))
		}
		for (j in seq_len(20)) {
			if (j %% 2 == 0) {
				maxi1 <- c(maxi1, range1[j])
			}
			if (j %% 2 != 0) {
				mini1 <- c(mini1, range1[j])
			}
		}
		for (j in seq_len(20)) {
			if (j %% 2 == 0) {
				maxi2 <- c(maxi2, range2[j])
			}
			if (j %% 2 != 0) {
				mini2 <- c(mini2, range2[j])
			}
		}
		for (r in seq_len(10)) {
			for (i in 1:(tamM + tamLNC)) {
				if (i <= tamM)
					datah[i, r] = (datah[i, r] - mini1[r]) / (maxi1[r] - mini1[r])
				if (i > tamM && i <= tamM + tamLNC)
					datah[i, r] = (datah[i, r] - mini2[r]) / (maxi2[r] - mini2[r])
			}
		}
	} else{
		for (u in seq_len(10)) {
			range1 <- c(range1, range(datah[1:tamM, u]))
		}
		for (u in seq_len(10)) {
			range2 <- c(range2, range(datah[(tamM + 1):(tamM + tamLNC), u]))
		}
		for (u in seq_len(10)) {
			range3 <-
				c(range3, range(datah[(tamM + tamLNC + 1):(tamM + tamLNC + tamSNC), u]))
		}
		for (j in seq_len(20)) {
			if (j %% 2 == 0) {
				maxi1 <- c(maxi1, range1[j])
			}
			if (j %% 2 != 0) {
				mini1 <- c(mini1, range1[j])
			}
		}
		for (j in seq_len(20)) {
			if (j %% 2 == 0) {
				maxi2 <- c(maxi2, range2[j])
			}
			if (j %% 2 != 0) {
				mini2 <- c(mini2, range2[j])
			}
		}
		for (j in seq_len(20)) {
			if (j %% 2 == 0) {
				maxi3 <- c(maxi3, range3[j])
			}
			if (j %% 2 != 0) {
				mini3 <- c(mini3, range3[j])
			}
		}
		for (r in seq_len(10)) {
			for (i in 1:(tamM + tamLNC + tamSNC)) {
				if (i <= tamM)
					datah[i, r] = (datah[i, r] - mini1[r]) / (maxi1[r] - mini1[r])
				if (i > tamM && i <= tamM + tamLNC)
					datah[i, r] = (datah[i, r] - mini2[r]) / (maxi2[r] - mini2[r])
				if (i > tamM + tamLNC)
					datah[i, r] = (datah[i, r] - mini3[r]) / (maxi3[r] - mini3[r])
			}
		}
	}
	return(datah)
}
