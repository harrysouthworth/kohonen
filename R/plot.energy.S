"plot.energy" <- 
function(object)
{
	if(length(object$energy$extend.ss) != 1) {
		breaks <- which(diff(object$energy$rayon) != 0)
		step <- object$energy$iteration
		extended.ss <- object$energy$extend.ss
		neighborhood <- c(0, breaks)
		plot(c(0, max(step)), range(extended.ss), type = "n", xlab = 
			"", ylab = "")
		lines(step, extended.ss)
		abline(v = neighborhood)
		axis(side = 3, tck = 0, at = neighborhood + diff(c(
			neighborhood, object$energy$lastiter))/2, labels = 
			paste(object$energy$rayon[neighborhood + 1], 
			"neighbours"), cex = 0.5)
		title(main = "Extended intra-inertia during learning process")
		title(xlab = "Iteration")
		title(ylab = "Extended intra-inertia")
	}
	else stop("no energy to plot")
	invisible()
}

