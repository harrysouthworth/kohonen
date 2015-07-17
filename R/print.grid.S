"print.grid" <- 
function(object)
{
	cat(" Grid layout for", object$nclass, "classes :", object$xdim, "*",
		object$ydim, "\n ")
	print(object$data)
}

