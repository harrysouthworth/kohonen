"get.grid.data" <- 
function(object, counts)
{
	computed.grid <- object$grid$coordinates
	if(any(is.na(counts))) {
		counts[is.na(counts)] <- 0
	}
	if(any(is.infinite(counts))) {
		counts[is.infinite(counts)] <- 0
	}
	if(object$grid$type == "hexagonal") {
		i <- as.hexagons(computed.grid, counts)
	}
	else if(object$grid$type == "rectangular") {
		i <- as.squares(computed.grid, counts)
	}
	else stop(paste("unimplemented type:", object$grid$type))
	i
}

