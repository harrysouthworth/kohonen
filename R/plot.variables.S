"plot.variables" <- 
function(object, name, ...)
{
	x <- object$weights
	dimnames(x)[[2]] <- names(object$data)
	noms <- c(names(object$data), names(object$list.contr))
	if(!all(is.element(name, noms)))
		stop(paste(name, "should be one of these values:", paste(
			noms, collapse = " ; ")))
	if(all(is.element(name, names(object$list.contr)))) {
		index <- which(is.element(name, names(object$list.contr)))
		name <- paste(names(object$list.contr)[index], dimnames(object$
			list.contr[[index]])[[1]], sep = "")
	}
	ylim <- range(x[, name])
	for(i in name)
		boxplot(x[, i], main = i, ylim = ylim, ...)
	invisible()
}

