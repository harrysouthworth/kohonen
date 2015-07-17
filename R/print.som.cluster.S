"print.som.cluster" <- 
function(object)
{
	cat("data subset:\n")
	if(object$is.empty)
		cat("\t<empty cluster>")
	else print(object$data)
	cat("\nneuron's values:\n")
	print(object$weights)
	invisible()
}

