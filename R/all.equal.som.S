"all.equal.som" <- 
function(target, current, tolerance = std.tolerance(), scale, ...)
{
	msg <- NULL
	lt <- length(target)
	lc <- length(current)
	if(lt != lc)
		return(paste("Numeric: lengths (", lt, ", ", lc, ") differ"),
			sep = "")
	else msg <- NULL
	target <- as.vector(target)
	current <- as.vector(current)
	out <- is.na(target)
	if(any(out != is.na(current)))
		return(paste("Missing value mismatches:", sum(is.na(current)),
			"in current,", sum(out), " in target"))
	out <- out | (target == current)
	if(all(out))
		return(T)
	target <- target[!out]
	current <- current[!out]
	xy <- mean(abs(target - current))
	if(missing(scale)) {
		xn <- mean(abs(target))
		if(xn > tolerance) {
			xy <- xy/xn
			what <- "relative"
		}
		else what <- "absolute"
	}
	else {
		xy <- xy/scale
		what <- "scaled"
	}
	if(is.na(xy) || xy > tolerance)
		F
	else T
}

