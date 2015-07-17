"boxplot.summary.som.cluster" <- 
function(object, plot = T, k = 1.5, boxWidth = 0.5)
{
	if(boxWidth > 1)
		stop("boxWidth must be <= 1")
	x <- object$ss
	if(object$is.empty && length(x) > 5)
		return(invisible())
	Quartils <- as.vector(summary(x)[, c("1st Qu.", "Median", "3rd Qu.")])
	inf.limit <- (x - Quartils[3])/(Quartils[3] - Quartils[1])
	sup.limit <- (Quartils[1] - x)/(Quartils[3] - Quartils[1])
	boolOut <- (inf.limit > k) | (sup.limit > k)
	ans <- which(boolOut)
	if(plot) {
		plot(c(-1, 1), range(x), type = "n", xaxt = "n", xlab = "",
			ylab = "")
		polygon(c( - boxWidth,  - boxWidth, boxWidth, boxWidth), c(
			Quartils[1], Quartils[3], Quartils[3], Quartils[1]),
			col = 6)
		lines(c( - boxWidth, boxWidth), rep(range(x[!boolOut])[1],
			2), col = 6, lwd = 2)
		lines(c( - boxWidth, boxWidth), rep(range(x[!boolOut])[2],
			2), col = 6, lwd = 2)
		lines(c( - boxWidth, boxWidth), rep(Quartils[2], 2), col = 0,
			lwd = 2)
		text(x = rep(0, sum(boolOut)), y = x[boolOut], labels = 
			as.character(dimnames(object$data[ans,  ])[[1]]),
			col = 8)
		lines(c(0, 0), c(Quartils[1], range(x[!boolOut])[1]), col = 6,
			lwd = 2, lty = 2)
		lines(c(0, 0), c(Quartils[3], range(x[!boolOut])[2]), col = 6,
			lwd = 2, lty = 2)
		polygon(c( - boxWidth,  - boxWidth, boxWidth, boxWidth), c(
			Quartils[1], Quartils[3], Quartils[3], Quartils[1]),
			density = 0, col = 6)
		key(text = list(c("Outliers", "Mean")), points = list(pch = 
			list(18, 16)), col = c(8, 5))
	}
	invisible(list(criteria = c((k * (Quartils[3] - Quartils[1])) + 
		Quartils[3], Quartils[1] - (k * (Quartils[3] - Quartils[1]))),
		data = object$data[ans,  , drop = F]))
}

