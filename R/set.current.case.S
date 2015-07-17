"set.current.case" <- 
function(object, current.index)
{
	if(current.index > object$grid$nclass || current.index < 1)
		stop("this index doesn't exist")
	object$grid$current.index <- current.index
	object
}

