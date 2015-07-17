"plot.som" <- 
function(object, type = "effectif", cex = 0.75, nom.variable = character(
	), precision = 4, cex.label = 0.75, ...)
{
	if(type == "variable") {
		if(length(nom.variable) == 0) {
			nom.variable <- c(names(object$data), names(object$
				list.contr))
			if(length(object$list.contr) > 0) {
				fact <- unlist(lapply(1:length(object$
					list.contr), function(x, object)
				{
					paste(names(object$list.contr)[x],
						dimnames(object$list.contr[[
						x]])[[1]], sep = "")
				}
				, object))
				nom.variable <- setdiff(nom.variable, fact)
			}
		}
		for(i in nom.variable)
			plot.variables(object, name = i, ...)
		return(invisible())
	}
	counts <- switch(type,
		distance = dist2full(dist(object$weights))[object$grid$
			current.index,  ],
		mean.ss = {
			object$withinss[object$size != 0] <- object$withinss[
				object$size != 0]/object$size[object$size !=
				0]
			object$withinss
		}
		,
		effectif = object$size,
		meta = object$metagroups,
		energy = {
			plot.energy(object)
			return(invisible())
		}
		,
		dxdy = {
			plot.dxdy(object, ...)
			return(invisible())
		}
		,
		stop("plot.som : unimplemented type"))
	if(type == "mean.ss") {
		counts[object$size == 0] <- 0
	}
	if(type == "variable") {
		counts[counts == NA] <- 0
	}
	counts <- round(counts, 8)
	i <- get.grid.data(object, counts)
	labels <- paste("class", 1:object$grid$nclass, "\n", signif(counts,
		precision))
	if(type != "meta")
		plot(i, labels, cex = cex, cex.label = cex.label)
	else plot(i, labels, cex = cex, cex.label = cex.label, col2use = sort(
			unique(counts)) + 1)
	cex <- cex * 1.5
	switch(type,
		distance = title(main = paste("Distances between class", 
			object$grid$current.index, "and each other"), sub = 
			"euclidean distances", cex = cex),
		effectif = title(main = paste("Size of the population"), sub
			 = "number of cases / class", cex = cex),
		meta = title(main = "Meta groups", cex = cex),
		mean.ss = title(main = "Means sums of squares", cex = cex),

			NULL)
	invisible()
}

