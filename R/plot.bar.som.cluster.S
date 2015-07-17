"plot.bar.som.cluster" <- 
function(object, col, alpha.x)
{
	plot(c(1, length(object)), c(0, 1), type = "n", xaxt = "n", yaxt = "n",
		ylab = "", xlab = "", bty = "n", axes = F)
	polygon(c(1, 1:length(object), length(object), 1), c(0, object, 0,
		NA), col = col, lwd = 1)
	invisible()
}

