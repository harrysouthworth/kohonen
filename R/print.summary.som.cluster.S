"print.summary.som.cluster" <- 
function(object)
{
	print.som.cluster(object)
	cat("\nMean sum of squares:\n")
	print(object$withinss)
	cat("\nStandart deviation sum of squares:\n")
	print(object$stdevss)
	cat("\nsums of squares:\n")
	print(object$ss)
	if(!object$is.empty && length(object$ss) > 5) {
		if(nrow(object$outliers) > 0) {
			cat("\nOutliers (method boxplot):\n")
			print(object$outliers)
		}
	}
	invisible()
}

