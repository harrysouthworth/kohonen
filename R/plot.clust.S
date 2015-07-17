"plot.clust" <- 
function(object, k, h, interactive = F)
{
	if(missing(k) && missing(h) && interactive == F)
		stop("choose a method for cuting the tree")
	ht <- hclust(dist(object$weights))
	plclust(tree = ht)
	if(interactive) {
		pointeur <- locator(1)
		meta.groups <- cutree(ht, h = pointeur$y)
	}
	else {
		pointeur <- list(x = 0, y = 0)
		if(missing(k) && !missing(h)) {
			meta.groups <- cutree(ht, h = h)
			pointeur$y <- h
		}
		if(!missing(k) && missing(h)) {
			meta.groups <- cutree(ht, k = k)
			pointeur$y <- max(attr(meta.groups, "height"))
		}
		if((!missing(k) && !missing(h)) || (missing(k) && missing(
			h)))
			stop("you must specify k or h")
	}
	segments(1, pointeur$y, nrow(object$weights), pointeur$y)
	object$metagroups <- as.vector(meta.groups)
	object
}

