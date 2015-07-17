"update.som" <- 
function(object, newdata)
{
	if(missing(newdata))
		dataset <- object$data
	else {
		dataset <- newdata
		if(length(object$list.contr) > 0) {
			dataset <- model.matrix(object$terms, dataset, object$
				list.contr)
			attributes(dataset) <- attributes(dataset)[c("dim",
				"dimnames")]
			dataset <- as.data.frame(dataset)
		}
	}
	if(!all(is.element(names(object$data), names(dataset))))
		stop("missing columns for updating the object")
	if(!all(is.element(names(dataset), names(object$data))))
		stop("too many columns in newdata")
	object$data <- dataset
	object$groups <- predict(object)[, 1]
	object$size <- size(object)
	withinss <- split(ss.per.row(object), object$groups)
	withinss <- sapply(withinss, function(x)
	sum(x))
	temp.withinss <- integer(object$grid$nclass)
	temp.withinss[as.integer(names(withinss))] <- withinss/object$size[
		as.integer(names(withinss))]
	object$withinss <- temp.withinss
	object
}

