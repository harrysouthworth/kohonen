"as.hexagons" <- 
function(coord, count)
{
	coord <- as.matrix(data.frame(x = coord[, 1], y = coord[, 2]))
	hex.prototype <- hexagone(rayon = 0.5, offset.angle = pi/6)
	dimensions <- lapply(hex.prototype, function(x)
	diff(range(x)))
	hex.coord <- .C("get_coordinates",
		as.vector(t(as.matrix(coord))),
		nrow(coord),
		hex.prototype$x,
		hex.prototype$y,
		double(nrow(coord) * 2 * 8),
		as.integer(7))
	hex.coord <- hex.coord[[5]]
	hex.coord[hex.coord == -1] <- NA
	xlim <- c(min(coord[, 1]) - dimensions[[1]]/2, max(c(coord[, 1] + 
		dimensions[[1]]/2, coord[, 2] + dimensions[[2]]/2)))
	ylim <- c(min(coord[, 2]) - dimensions[[2]]/2, max(c(coord[, 1] + 
		dimensions[[1]]/2, coord[, 2] + dimensions[[2]]/2)))
	ans <- list(polygones = matrix(hex.coord, ncol = 2, byrow = T), count
		 = count, xlim = xlim, ylim = ylim, coord = coord, proto = 
		hex.prototype)
	class(ans) <- "map"
	ans
}

