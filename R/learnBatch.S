"learnBatch" <- 
function(object, number.iter, max.rayon = 3, min.rayon = 0)
{
	nrow.dataset <- dim(object$data)[1]
	ncol.dataset <- dim(object$data)[2]
	if(missing(number.iter))
		stop("parameter 'number.iter' must be defined")
	if(max.rayon < min.rayon)
		stop("max.rayon must be > min.rayon")
	if(min.rayon < 0)
		stop("min.rayon must be >= 0")
	rayon <- as.integer(ceiling(seq(max.rayon, min.rayon - 
		0.98999999999999999, length = number.iter)))
	vector.dataset <- as.vector(t(as.matrix(object$data)))
	vector.weigths <- as.vector(t(object$weights))
	computed.dist.coordinates <- as.vector(t(object$grid$
		distances.coordinates))
	dist.step <- matrix(object$grid$distances.coordinates[1,  ], byrow = T,
		ncol = object$grid$ydim)
	if(object$grid$xdim == 1)
		dist.step <- dist.step[1, 2]
	else if(object$grid$ydim == 1)
		dist.step <- dist.step[2, 1]
	else if(object$grid$type == "hexagonal") {
		dist.step <- dist.step[1, 2]
	}
	else {
		dist.step <- diag(dist.step[c(1:min(object$grid$ydim, object$
			grid$xdim)), c(1:min(object$grid$ydim, object$grid$
			xdim))])[2]
	}
	answer <- .C("KohonenCBatch",
		vector.dataset,
		nrow.dataset,
		ncol.dataset,
		vector.weigths,
		object$grid$nclass,
		rayon,
		as.integer(number.iter),
		computed.dist.coordinates,
		dist.step,
		integer(object$grid$nclass),
		integer(nrow.dataset),
		ifelse(object$neighborhood == "gaussian", 1, 0),
		double(object$grid$nclass),
		double(object$grid$nclass),
		double(number.iter))
	object$weights <- matrix(answer[[4]], byrow = T, ncol = ncol.dataset)
	object$energy$extend.ss <- append(object$energy$extend.ss, answer[[
		15]])
	object$energy$rayon <- append(object$energy$rayon, rayon)
	object$energy$iteration <- append(object$energy$iteration, (object$
		energy$lastiter + 1):(object$energy$lastiter + number.iter))
	object$energy$lastiter <- number.iter + object$energy$lastiter
	object <- update.som(object)
	object$n.iterations <- number.iter
	object
}

