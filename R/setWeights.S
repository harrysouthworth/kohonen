"setWeights" <- 
function(x, new)
{
	if(class(x) != "som")
		stop("class of x must be 'som'")
	if(class(new) != "matrix")
		stop("class of x must be 'matrix'")
	if(!all(dim(x$weights) == dim(new)))
		stop("class of x must be 'matrix'")
	x$weights <- new
	x <- update(x)
	x
}

