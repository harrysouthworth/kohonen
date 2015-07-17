"grid" <- 
function(xdim, ydim, type = "hexagonal", current.index = 1)
{
	val <- NULL
	if(length(xdim) != 1 || length(ydim) != 1 || !is.integer(xdim) || !
		is.integer(ydim) || xdim <= 0 || ydim <= 0)
		stop("xdim and ydim must be positive non null integers of length 1"
			)
	if(length(type) != 1 && !is.element(type, c("hexagonal", "rectangular"
		)))
		stop("type must be 'hexagonal' or 'rectangular'")
	if(current.index <= 0 || current.index > (xdim * ydim))
		stop(paste("current.index must be >= 1 and <=", (xdim * ydim))
			)
	grid.val <- matrix(1:(xdim * ydim), nrow = xdim, byrow = T)
	if(ydim == xdim && ydim == 1)
		stop("xdim and ydim must be greater than 1")
	if(ydim == 1)
		grid.val <- matrix(grid.val)
	else if(xdim == 1)
		grid.val <- matrix(grid.val, nrow = 1)
	else grid.val <- grid.val[xdim:1,  ]
	computed.grid <- compute.grid(grid.val, type)
	distances.computed.grid <- dist2full(dist(computed.grid))
	val <- list("grid", data = grid.val, type = type, xdim = xdim, ydim = 
		ydim, nclass = xdim * ydim, coordinates = computed.grid, 
		distances.coordinates = distances.computed.grid, counts = seq(
		1:(xdim * ydim)), current.index = current.index)
	class(val) = "grid"
	val
}

