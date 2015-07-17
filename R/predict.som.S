"predict.som" <- 
function(object, newdata, type = "class")
{
	#--- predict data based on model and dataset
	if(missing(newdata)) dataset <- object$data else {
		dataset <- newdata
		if(length(object$list.contr) > 0) {
			dataset <- model.matrix(object$terms, dataset, object$
				list.contr)
		}
		attributes(dataset) <- attributes(dataset)[c("dim", "dimnames"
			)]
		dataset <- as.data.frame(dataset)
	}
	dist.per.class <- .C("rows_dist_weights",
		as.vector(t(dataset)),
		nrow(dataset),
		ncol(dataset),
		as.vector(t(object$weights)),
		object$grid$nclass,
		double(nrow(dataset) * object$grid$nclass))[[6]]
	dist.per.class <- matrix(dist.per.class, nrow = nrow(dataset), byrow
		 = T)
	dimnames(dist.per.class)[[2]] <- paste("cluster", 1:object$grid$nclass,
		sep = "")
	output <- list()
	for(i in type) {
		if(i == "class" || i == "meta") {
			answer <- apply(dist.per.class, 1, FUN = function(x)
			which(x == min(x)))
			if(i == "meta") {
				answer <- object$metagroups[answer]
			}
			output[[i]] <- answer
		}
		else if(i == "dist") {
			output[[i]] <- dist.per.class
		}
		else stop("type must be class, meta or dist")
	}
	answer <- as.data.frame(output)
	row.names(answer) <- row.names(dataset)
	answer
}

