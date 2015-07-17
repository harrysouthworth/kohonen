"print.summary.som" <- 
function(object)
{
	cat("inertia:\n")
	print(object$inertia)
	cat("\nClustering vector:\n")
	print(as.vector(object$class[, 1]))
	cat("\nMeta clustering vector:\n")
	print(as.vector(object$metaclass[, 1]))
	cat("\nWithin cluster sum of squares:\n")
	print(object$withinss)
	cat("\nCluster sizes:\n")
	print(object$size)
	cat("\nAvailable arguments:\n")
	print(names(object))
}

