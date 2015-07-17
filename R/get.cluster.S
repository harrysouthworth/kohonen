"get.cluster" <- 
function(object, class.number)
{
	ans <- list(data = NULL, weights = NULL, class.number = class.number,
		withinss = NULL, ss = NULL, is.empty = T)
	dataset <- clean.dataset(object$data)
	weights <- clean.dataset(object$weights)
	clusters <- object$groups
	clusters <- which(clusters == ans$class.number)
	ans$weights <- as.vector(weights[class.number,  ])
	names(ans$weights) <- names(dataset)
	temp <- as.matrix(rbind(dataset, weights))
	.range <- apply(temp, 2, range)
	temp <- reduce2021(temp, range = .range)
	ans$ylim <- range(temp)
	ans$weights2plot <- clean.dataset(reduce2021(weights[class.number,
		, drop = F], range = .range))
	temp <- clean.dataset(as.matrix(rbind(ans$data2plot, ans$weights2plot)
		))
	if(length(clusters) > 0) {
		ans$is.empty <- F
		ans$data <- object$data[clusters,  , drop = F]
		withinss <- ss.per.row(object)[clusters]
		ans$ss <- withinss
		ans$withinss <- object$withinss[ans$class.number]
		ans$data2plot <- clean.dataset(reduce2021(ans$data, range = 
			.range))
		dimnames(ans$weights2plot)[[2]] <- names(dataset)
	}
	oldClass(ans) <- "som.cluster"
	ans
}

