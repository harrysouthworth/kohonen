"print.som" <- 
function(object)
{
	cat("------som object------\n")
	cat("Grid layout for", object$grid$nclass, "classes :", object$grid$
		xdim, "*", object$grid$ydim, "\n")
	cat("Type of grid is :", object$grid$type, "\n")
	if(object$n.iterations < 1) {
		cat("This object has never been trained. Use the function 'learn' or 'learnBatch'.\n"
			)
	}
	else {
		cat("Mean of within sums of squares :", signif(mean(object$
			withinss), 5), "\n")
		cat("Class size :\n")
		print(object$size)
	}
	invisible()
}

