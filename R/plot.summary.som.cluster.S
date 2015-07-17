"plot.summary.som.cluster" <- 
function(object, no.annotation = F)
{
	if(no.annotation)
		plot(c(1, ncol(object$weights2plot)), object$ylim, type = "n",
			xaxt = "n", ylab = "", xlab = "", bty = "n", axes = F)
	else plot(c(1, ncol(object$weights2plot)), c(0, 1), type = "n", xaxt
			 = "n", ylab = "", xlab = "", main = paste(
			"data values for class", object$class.number))
	key.text <- c(paste("data (N = ", nrow(object$data2plot), ")", sep = 
		""), "mean", "median", "prototype")
	if(!is.null(object$data2plot)) {
		try(apply(object$data2plot, 1, FUN = function(x)
		matlines(1:length(x), x, lwd = 1)))
	}
	matlines(1:ncol(object$weights2plot), as.vector(object$weights2plot),
		lty = 3, lwd = 2)
	if(!no.annotation) {
		key(text = list(c(paste("data (N = ", nrow(object$data2plot),
			")", sep = ""), "prototype"), adj = 1), line = list(
			lty = c(1, 3)), transparent = T)
	}
	invisible()
}

