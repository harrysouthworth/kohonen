"as.squares" <- 
function(coord, count)
{
	coord <- as.matrix(data.frame(x = coord[, 1], y = coord[, 2]))
	prototype <- square(rayon = 0.5)
	dimensions <- list(x = 1, y = 1)
	hex.coord <- .C("get_coordinates",
		as.double(as.vector(t(as.matrix(coord)))),
		nrow(coord),
		prototype$x,
		prototype$y,
		double(nrow(coord) * 2 * 5),
		as.integer(4))
	hex.coord <- hex.coord[[5]]
	hex.coord[hex.coord == -1] <- NA
	xlim <- c(min(coord[, 1]) - dimensions[[1]]/2, max(c(coord[, 1] + 
		dimensions[[1]]/2, coord[, 2] + dimensions[[2]]/2)))
	ylim = c(min(coord[, 2]) - dimensions[[2]]/2, max(c(coord[, 1] + 
		dimensions[[1]]/2, coord[, 2] + dimensions[[2]]/2)))
	ans <- list(polygones = matrix(hex.coord, ncol = 2, byrow = T), count
		 = count, xlim = xlim, ylim = ylim, coord = coord, proto = 
		prototype)
	class(ans) <- "map"
	ans
}

