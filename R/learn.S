"learn" <- 
function(object, number.iter, order, max.alpha = 0.050000000000000003, 
	min.alpha = 0., max.rayon = 3, min.rayon = 0, step.eval.si = 100)
{
	nrow.dataset <- dim(object$data)[1]
	ncol.dataset <- dim(object$data)[2]
	if(missing(number.iter) && missing(order))
		stop("parameter 'order' or 'number.iter' must be defined")
	if(!missing(number.iter))
		order <- rsample(n = nrow.dataset, replace = T, size = 
			number.iter)
	else number.iter <- length(order)
	if(!all(is.element(order, 1:nrow.dataset)))
		stop("order contains wrong index for the defined dataset")
	if(max.rayon < min.rayon)
		stop("max.rayon must be > min.rayon")
	if(min.rayon < 0)
		stop("min.rayon must be >= 0")
	if(max.alpha < min.alpha)
		stop("max.alpha must be > min.alpha")
	if(max.alpha > 1)
		stop("max.alpha must be <= 1")
	if(min.alpha < 0)
		stop("min.alpha must be >= 0")
	rayon <- as.integer(ceiling(seq(max.rayon, min.rayon - 
		0.98999999999999999, length = number.iter)))
	alpha <- seq(max.alpha, min.alpha, length = number.iter)
	vector.dataset <- as.vector(t(as.matrix(object$data)))
	vector.weigths <- as.vector(t(object$weights))
	computed.dist.coordinates <- as.vector(t(object$grid$
		distances.coordinates))
	order <- sample(1:nrow.dataset, replace = T, size = number.iter)
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
	if(step.eval.si >= 1)
		step.eval.si <- as.integer(c(seq(0, number.iter - 1, by = 
			step.eval.si), number.iter))
	else step.eval.si <- integer(1)
	answer <- .C("KohonenC",
		vector.dataset,
		nrow.dataset,
		ncol.dataset,
		vector.weigths,
		object$grid$nclass,
		alpha,
		rayon,
		length(rayon),
		computed.dist.coordinates,
		order,
		dist.step,
		ifelse(object$neighborhood == "gaussian", 1, 0),
		double(object$grid$nclass),
		integer(nrow.dataset),
		double(object$grid$nclass),
		step.eval.si,
		double(length(step.eval.si)))
	object$weights <- matrix(answer[[4]], byrow = T, ncol = ncol.dataset)
	object$energy$extend.ss <- append(object$energy$extend.ss, answer[[
		17]])
	object$energy$rayon <- append(object$energy$rayon, rayon)
	object$energy$iteration <- append(object$energy$iteration, 
		step.eval.si + object$energy$lastiter)
	object$energy$lastiter <- number.iter + object$energy$lastiter
	object$range.rayon <- c(min.rayon, max.rayon)
	object$range.alpha <- c(min.alpha, max.alpha)
	object <- update.som(object)
	object$n.iterations <- number.iter
	object
}

