"summary.som.cluster" <- 
function(object)
{
	Names <- dimnames(object$data)[[1]]
	object$outliers <- NULL
	object$stdevss <- NULL
	if(!object$is.empty) {
		if(length(Names) > 5)
			object$outliers <- boxplot.summary.som.cluster(object,
				plot = F, k = 1.5)$data
		object$stdevss <- stdev(object$ss)
	}
	oldClass(object) <- "summary.som.cluster"
	object
}

