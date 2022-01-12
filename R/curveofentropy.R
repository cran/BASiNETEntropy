#'@title Creates an entropy curve
#'@name curveofentropy
#'
#'@description A function that from the entropy measures and threshold creates
#'an entropy curve.
#'@param H The 'training' return for the entropy measures
#'@param threshold The 'training' return for the threshold
#'
#'@return Returns a entropy curve
#'@author Murilo Montanini Breve
#'@importFrom graphics abline axis legend text plot
curveofentropy <- function(H, threshold) {
	H <- unlist(H)
	threshold <- unlist(threshold)
	plot(H,
	     main = "Curve of entropy",
	     ylab = "Sum of Entropies",
	     xlab = "Edges distribution",
	     xaxt = "n")
	axis(1, xaxp = c(0, 4096, 8), las = 1)
	abline(v = threshold, col = "red")
	text(print(paste0("T = ", threshold)),
	     x = threshold,
	     y = 0,
	     srt = 90,
	     adj = c(0,-0.5))
	legend("topright",
	       legend = "Threshold",
	       pch = "|",
	       col = "red")
}
