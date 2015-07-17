"compute.grid" <- 
function(object, type)
{
	if(type == "rectangular") {
		y <- 1:dim(object)[1]
		x <- 1:dim(object)[2]
		coord <- as.matrix(expand.grid(x = x, y = y))
	}
	else {
		hex.prototype <- hexagone(rayon = 0.5, offset.angle = pi/
			6)
		dimensions <- lapply(hex.prototype, function(x)
		diff(range(x)))
		list.coord <- lapply(1:nrow(object), function(i, dims, sizes)
		{
			ans <- list()
			if(i %% 2 == 1)
				ans[["x"]] <- seq(sizes[[1]]/2, sizes[[1]] *
					dims[2] - sizes[[1]]/2, by = sizes[[
					1]])
			else ans[["x"]] <- seq(sizes[[1]], sizes[[1]] * dims[
					2] + sizes[[1]]/2, by = sizes[[1]])
			ans[["y"]] <- sizes[[2]]/2 + (i - 1) * (sizes[[2]] -
				0.25)
			ans
		}
		, dim(object), dimensions)
		coord <- as.matrix(data.frame(x = unlist(lapply(list.coord,
			function(x)
		x$x)), y = unlist(lapply(list.coord, function(x, n)
		rep(x$y, n), dim(object)[2]))))
	}
	coord
}

