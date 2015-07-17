"summary.som" <- 
function(object)
{
	answer <- list()
	answer$inertia <- c(inter.inertia(object), intra.inertia(object))
	answer$weights <- object$weights
	answer$withinss <- object$withinss
	answer$size <- object$size
	answer$class <- predict(object)
	answer$metaclass <- predict(object, type = "meta")
	names(answer$inertia) <- c("inter", "intra")
	class(answer) <- "summary.som"
	answer
}

