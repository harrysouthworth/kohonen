"plot.map" <- 
function(i, labels, cex = 1, cex.label, col2use)
{
	distances <- i$count
	if(missing(col2use)) {
		minx <- min(distances)
		maxx <- max(distances)
		distances <- cut(distances, seq(minx, maxx, length = 100),
			include.lowest = T)
		if(use.device.palette()) {
			tp <- trellis.par.get("regions")$col
		} else {
			tp <- image.palette()
		}
	}
	else {
		tp <- col2use
	}
	ratio.x <- 0.85
	ratio.y <- 1
	i$polygones[, 1] <- i$polygones[, 1] * ratio.x
	i$polygones[, 2] <- i$polygones[, 2] * ratio.y
	i$coord[, "x"] <- i$coord[, "x"] * ratio.x
	i$coord[, "y"] <- i$coord[, "y"] * ratio.y
	i$proto$x <- i$proto$x * ratio.x
	i$proto$y <- i$proto$y * ratio.y
	plot(bty = "n", xlab = "", ylab = "", axes = F, i$xlim, i$ylim, type = 
		"n")
	polygon(i$polygones[, 1], i$polygones[, 2], col = tp[distances], border
		 = F)
	lines(i$polygones[, 1], i$polygones[, 2])
	if(cex.label > 0.001)
		text(i$coord[, 1], i$coord[, 2], labels, cex = cex.label)
	n <- dim(i$coord)[1]
	if(all.equal.som(round(i$count, 0), i$count) && !missing(col2use)) {
		nint <- length(unique(i$count))
		if(length(col2use) > 1)
			key(x = i$xlim[2] * (ratio.x + 0.10000000000000001),
				y = i$ylim[2], text.width.multiplier = 1, 
				corner = c(0, 1), transparent = T, rectangle = 
				list(size = 2, col = tp[(1:nint)]), text = list(
				as.character(sort(unique(i$count))), col = 1,
				adj = 0), between = 0.5, cex = cex)
	}
	else {
		nint <- 50
		x <- c(i$xlim[2] * (ratio.x + 0.10000000000000001), i$xlim[
			2])
		y <- c(i$ylim[1], i$ylim[2])
		int <- interp(i$coord[, 1] + rnorm(n, 0, 0.001), i$coord[,
			2] + rnorm(n, 0, 0.001), as.double(i$count))
		image.legend(int, horizontal = F, x = x, y = y, nint = nint,
			cex = cex)
	}
	invisible(i)
}

