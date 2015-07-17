"inter.inertia" <- 
function(object)
{
	answer <- .C("get_inter_ss",
		as.vector(t(object$weights)),
		nrow(object$weights),
		ncol(object$weights),
		as.vector(apply(object$data, 2, median)),
		double(object$grid$nclass))[[5]]
	sum(answer * object$size)
}

