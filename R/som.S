"som" <- 
function(formula, data, grid = grid(9, 9, type = "hexagonal"), neighborhood = 
	"gaussian", weights.min = -1, weights.max = 1)
{
	m <- match.call(expand.dots = F)
	call <- match.call()
	m$grid <- m$neighborhood <- m$x <- m$weights.matrix <- m$contrasts <-
		m$weights.min <- m$weights.max <- NULL
	m$drop.unused.levels <- F
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	Terms <- attr(m, "terms")
	xvars <- as.character(attr(Terms, "variables"))
	Terms <- attr(m, "terms")
	attr(Terms, "intercept") <- 0
	if((yvar <- attr(Terms, "response")) > 0)
		stop("unsupervised model can not have any response. see formula argument."
			)
	if(length(xvars) > 0) {
		xlevels <- lapply(m[xvars], levels)
		xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
		if(length(xlevels) == 0)
			xlevels <- NULL
	}
	else xlevels <- NULL
	m.contrasts <- list()
	contrasts <- list()
	wf <- which(sapply(m, is.factor))
	if(length(wf) > 0) {
		contrasts <- as.list(rep("contr.treatment", length(wf)))
		names(contrasts) <- names(m)[wf]
	}
	if(length(contrasts) > 0) {
		for(i in 1:length(contrasts)) {
			m.contrasts[[names(contrasts)[i]]] <- do.call(
				contrasts[[i]], list(n = levels(data[, names(
				contrasts)[i]])))
		}
		data <- model.matrix(Terms, m, m.contrasts)
	}
	else data <- model.matrix(Terms, m, data)
	attributes(data) <- attributes(data)[c("dim", "dimnames")]
	class(data) <- "matrix"
	m <- list(data = as.data.frame(data))
	weights.matrix <- matrix(runif(grid$nclass * ncol(data), min = 
		weights.min, max = weights.max), nrow = grid$nclass)
	var.names <- dimnames(data)[[2]]
	dimnames(weights.matrix) <- list(NULL, var.names)
	m$withinss <- rep(-1, grid$nclass)
	m$weights <- weights.matrix
	m$grid <- grid
	m$size = rep(-1, grid$nclass)
	m$neighborhood = switch(neighborhood,
		gaussian = "gaussian",
		uniform = "uniform",
		stop(paste("neighborhood", neighborhood, "not implemented")))
	m$groups = rep(-1, nrow(data))
	m$metagroups = rep(1, grid$nclass)
	m$terms <- Terms
	m$call <- call
	m$list.contr <- m.contrasts
	m$n.iterations <- 0
	m$range.rayon <- c(0, 0)
	m$range.alpha <- c(0, 0)
	m$energy <- list(lastiter = 0, iteration = integer(0), extend.ss = 
		numeric(0), rayon = integer(0))
	oldClass(m) <- "som"
	m
}

