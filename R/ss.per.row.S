"ss.per.row" <- 
function(object)
{
	answer <- .C("get_individual_ss",
		as.vector(t(object$data)),
		nrow(object$data),
		ncol(object$data),
		as.vector(t(object$weights)),
		nrow(object$weights),
		as.vector(object$groups),
		double(nrow(object$data)))
	as.vector(answer[[7]])
}

