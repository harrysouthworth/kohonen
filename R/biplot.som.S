"biplot.som" <- 
function(object, type = "lines", cex = 0.75)
{
	counts <- object$metagroups
	legend = T
	if(length(unique(counts)) == 1) {
		counts <- rep(1, length(object$metagroups))
		legend = F
	}
	spec.grid <- get.grid.data(object, counts)
	labels <- rep("", object$grid$nclass)
	spec.grid <- plot.map(spec.grid, labels, cex = cex, cex.label = 0,
		col2use = sort(unique(counts)) + 1)
	main <- paste(type, "biplot")
	if(object$grid$type == "hexagonal") {
		list.squ <- list(x = rep(c(spec.grid$proto$x[c(1, 3, 4, 6)],
			NA), object$grid$nclass) + rep(spec.grid$coord[, "x"],
			each = 5), y = rep(c(spec.grid$proto$y[c(1, 3, 4, 6)],
			NA), object$grid$nclass) + rep(spec.grid$coord[, "y"],
			each = 5))
		list.mat.squ <- list(x = matrix(rep(spec.grid$proto$x[c(1,
			4)], object$grid$nclass) + rep(spec.grid$coord[, "x"],
			each = 2), ncol = 2, byrow = T), y = matrix(rep(
			spec.grid$proto$y[c(1, 4)], object$grid$nclass) + rep(
			spec.grid$coord[, "y"], each = 2), ncol = 2, byrow = T)
			)
	}
	else {
		list.squ <- list(x = rep(c(spec.grid$proto$x, NA), object$
			grid$nclass) + rep(spec.grid$coord[, "x"], each = 5),
			y = rep(c(spec.grid$proto$y, NA), object$grid$nclass) +
			rep(spec.grid$coord[, "y"], each = 5))
		list.mat.squ <- list(x = matrix(rep(spec.grid$proto$x[c(1,
			4)], object$grid$nclass) + rep(spec.grid$coord[, "x"],
			each = 2), ncol = 2, byrow = T), y = matrix(rep(
			spec.grid$proto$y[c(2, 4)], object$grid$nclass) + rep(
			spec.grid$coord[, "y"], each = 2), ncol = 2, byrow = T)
			)
	}
	if(type == "lines") {
		for(i in 1:object$grid$nclass) {
			subplot(x = list.mat.squ$x[i,  ], y = list.mat.squ$
				y[i,  ], fun = plot.summary.som.cluster(
				get.cluster(object, i), no.annotation = T))
		}
		key(spec.grid$xlim[1], text.width.multiplier = 1, corner = c(
			0, 0), transparent = T, text = list(c("data", 
			"prototype"), adj = 0), lines = list(lty = c(1, 3)),
			cex = cex, columns = 2)
	}
	else if(type == "solid") {
		object$weights <- reduce2021(object$weights, range = apply(
			object$weights, 2, range))
		for(i in 1:object$grid$nclass) {
			subplot(x = list.mat.squ$x[i,  ], y = list.mat.squ$
				y[i,  ], fun = plot.bar.som.cluster(object$
				weights[i,  ], col = trellis.settings$add.text$
				col))
		}
	}
	else if(type == "stars") {
		object$weights <- reduce2021(object$weights, range = apply(
			object$weights, 2, range))
		if(object$grid$type == "hexagonal")
			inch = max(spec.grid$proto$y[c(1, 4)]) * 2 * par("uin")[
				2]
		else {
			inch1 = max(spec.grid$proto$y[c(2, 4)]) * par("uin")[
				2]
			inch2 = max(spec.grid$proto$x[c(2, 4)]) * par("uin")[
				1]
			inch = min(inch1, inch2)
		}
		symbols(x = spec.grid$coord[, "x"], y = spec.grid$coord[, "y"],
			stars = object$weights, inches = inch, add = T, col = 
			trellis.settings$add.text$col)
	}
	else if(type == "boxplot") {
		ylim <- range(ss.per.row(object))
		for(i in 1:object$grid$nclass) {
			cl <- get.cluster(object, i)
			subplot(x = list.mat.squ$x[i,  ], y = list.mat.squ$
				y[i,  ], fun = boxplot(cl$ss, ylim = ylim,
				xaxt = "n", yaxt = "n", medlwd = 1, axes = F))
		}
	}
	title(main = main)
	invisible()
}

