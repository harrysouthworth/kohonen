"getWeights" <- 
function(x)
{
	if(class(x) != "som")
		stop("class of x must be 'som'")
	x$weights
}

