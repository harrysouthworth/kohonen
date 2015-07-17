"plot.dxdy" <- 
function(object, sample.size = 300, plot = T, col.points = 1, cex.points = 1,
	col.line = 1, pch.points = 1)
{
	fitted.vect <- fitted(object)
	id1 <- sample(1:length(fitted.vect), size = sample.size, replace = T)
	id2 <- sample(1:length(fitted.vect), size = sample.size, replace = T)
	bool.ok <- id1 != id2
	id1 <- id1[bool.ok]
	id2 <- id2[bool.ok]
	bool.ok <- duplicated(paste(id1, "-", id2, sep = ""))
	id1 <- id1[!bool.ok]
	id2 <- id2[!bool.ok]
	id.matrix <- matrix(c(id1, id2), byrow = F, ncol = 2)
	dimnames(id.matrix)[[2]] <- c("id1", "id2")
	class.matrix <- t(apply(id.matrix, 1, FUN = function(row, fitted)
	{
		fitted[row]
	}
	, fitted.vect))
	bool.ok <- class.matrix[, 1] != class.matrix[, 2]
	class.matrix <- class.matrix[bool.ok,  ]
	id.matrix <- id.matrix[bool.ok,  ]
	class.matrix <- class.matrix - 1
	id.matrix <- id.matrix - 1
	ans <- .C("arrangedist",
		id.matrix[, 1],
		id.matrix[, 2],
		class.matrix[, 1],
		class.matrix[, 2],
		as.vector(t(object$data)),
		as.vector(t(object$weights)),
		nrow(id.matrix),
		ncol(object$weights),
		double(2 * nrow(id.matrix)))
	DataToPlot <- matrix(ans[[9]], ncol = 2, byrow = T, dimnames = list(
		NULL, c("dx", "dy")))
	if(plot) {
		plot(list(x = range(DataToPlot[, 1]), y = range(DataToPlot[
			, 2])), type = "n", xlab = "n", ylab = "n")
		matpoints(DataToPlot[, 1], DataToPlot[, 2], type = "p", pch = 
			pch.points, col = col.points, cex = cex.points)
		fit.lm <- lm(DataToPlot[, 2] ~ DataToPlot[, 1])
		#abline( fit.lm , col = 2)
		abline(a = 0, b = 1)
		title("DY versus DX", xlab = "rows-between distances", ylab = 
			"class-between distances", axes = T)
	}
	invisible(cbind(DataToPlot, id.matrix))
}

