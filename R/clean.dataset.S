"clean.dataset" <- 
function(x)
{
	x[abs(x) < .Machine$double.eps] <- 0
	x[is.na(x)] <- 0
	x
}

