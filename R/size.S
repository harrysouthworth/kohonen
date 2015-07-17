"size" <- 
function(object)
{
	size <- rep(0, object$grid$nclass)
	names(size) <- 1:object$grid$nclass
	temp.size <- table(predict(object)[, 1])
	size[names(temp.size)] <- temp.size
	size
}

